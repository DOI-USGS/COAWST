      SUBROUTINE ad_uv3dmix4 (ng, tile)
!
!svn $Id: ad_uv3dmix4_geo.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes adjoint biharmonic mixing of momentum,        !
!  rotated along geopotentials,  from the horizontal divergence        !
!  of the stress tensor. A transverse isotropy is assumed so the       !
!  stress tensor is split into vertical and horizontal subtensors.     !
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
      USE mod_ncparam
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
      CALL wclock_on (ng, iADM, 33)
#endif
      CALL ad_uv3dmix4_tile (ng, tile,                                  &
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
     &                       GRID(ng) % ad_Hz,                          &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % ad_z_r,                         &
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
     &                       MIXING(ng) % Uvis3d_r,                     &
     &                       MIXING(ng) % Vvis3d_r,                     &
     &                       MIXING(ng) % ad_Uvis3d_r,                  &
     &                       MIXING(ng) % ad_Vvis3d_r,                  &
# else
     &                       MIXING(ng) % visc3d_r,                     &
     &                       MIXING(ng) % ad_visc3d_r,                  &
# endif
#else
     &                       MIXING(ng) % visc4_p,                      &
     &                       MIXING(ng) % visc4_r,                      &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                       DIAGS(ng) % DiaRUfrc,                      &
!!   &                       DIAGS(ng) % DiaRVfrc,                      &
!!   &                       DIAGS(ng) % DiaU3wrk,                      &
!!   &                       DIAGS(ng) % DiaV3wrk,                      &
#endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       COUPLING(ng) % ad_rufrc,                   &
     &                       COUPLING(ng) % ad_rvfrc)
#ifdef PROFILE
      CALL wclock_off (ng, iADM, 33)
#endif
      RETURN
      END SUBROUTINE ad_uv3dmix4
!
!***********************************************************************
      SUBROUTINE ad_uv3dmix4_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
#ifdef MASKING
     &                             pmask, rmask, umask, vmask,          &
#endif
     &                             om_p, om_r, om_u, om_v,              &
     &                             on_p, on_r, on_u, on_v,              &
     &                             pm, pn,                              &
     &                             Hz, ad_Hz,                           &
     &                             z_r, ad_z_r,                         &
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
     &                             Uvis3d_r, Vvis3d_r,                  &
     &                             ad_Uvis3d_r, ad_Vvis3d_r,            &
# else
     &                             visc3d_r, ad_visc3d_r,               &
# endif
#else
     &                             visc4_p, visc4_r,                    &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                             DiaRUfrc, DiaRVfrc,                  &
!!   &                             DiaU3wrk, DiaV3wrk,                  &
#endif
     &                             u, v, ad_u, ad_v,                    &
     &                             ad_rufrc, ad_rvfrc)
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
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
      real(r8), intent(in) :: Uvis3d_r(LBi:,LBj:,:)
      real(r8), intent(in) :: Vvis3d_r(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: visc3d_r(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: visc4_p(LBi:,LBj:)
      real(r8), intent(in) :: visc4_r(LBi:,LBj:)
# endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
      real(r8), intent(inout) :: ad_Uvis3d_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Vvis3d_r(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: ad_visc3d_r(LBi:,LBj:,:)
#  endif
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
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
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
      real(r8), intent(in) :: Uvis3d_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Vvis3d_r(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: visc3d_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
# else
      real(r8), intent(in) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc4_r(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
      real(r8), intent(inout) :: ad_Uvis3d_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Vvis3d_r(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: ad_visc3d_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, kk, kt, k1, k1b, k2, k2b

      real(r8) :: cff, fac1, fac2, pm_p, pn_p
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: cff5, cff6, cff7, cff8
      real(r8) :: dmUdz, dnUdz, dmVdz, dnVdz
#ifdef VISC_3DCOEF
      real(r8) :: Uvis_p, Vvis_p, visc_p
      real(r8) :: ad_fac1, ad_fac2,ad_Uvis_p, ad_Vvis_p, ad_visc_p
#endif
      real(r8) :: adfac, ad_cff
      real(r8) :: adfac1, adfac2, adfac3, adfac4, adfac5, adfac6
      real(r8) :: ad_cff1, ad_cff2, ad_cff3, ad_cff4
      real(r8) :: ad_cff5, ad_cff6, ad_cff7, ad_cff8
      real(r8) :: ad_dmUdz, ad_dnUdz, ad_dmVdz, ad_dnVdz

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: LapV

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ad_LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ad_LapV

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: UFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: UFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: VFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: VFsx
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

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_UFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_UFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_VFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_VFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dmUde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dmVde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dnUdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dnVdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dVdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZde_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZdx_r

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize private adjoint variables and arrays.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_cff4=0.0_r8
      ad_cff5=0.0_r8
      ad_cff6=0.0_r8
      ad_cff7=0.0_r8
      ad_cff8=0.0_r8

#ifdef VISC_3DCOEF
      ad_fac1=0.0_r8
      ad_fac2=0.0_r8
      ad_Uvis_p=0.0_r8
      ad_Vvis_p=0.0_r8
      ad_visc_p=0.0_r8
#endif

      ad_dmUdz=0.0_r8
      ad_dnUdz=0.0_r8
      ad_dmVdz=0.0_r8
      ad_dnVdz=0.0_r8

      ad_UFe(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_UFx(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_VFe(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_VFx(IminS:ImaxS,JminS:JmaxS)=0.0_r8

      ad_UFse(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_UFsx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_VFse(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_VFsx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8

      ad_dmUde(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dmVde(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dnUdx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dnVdx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8

      ad_dUdz(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dVdz(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8

      ad_dZde_p(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dZde_r(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dZdx_p(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dZdx_r(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8

      ad_LapU(IminS:ImaxS,JminS:JmaxS,1:N(ng))=0.0_r8
      ad_LapV(IminS:ImaxS,JminS:JmaxS,1:N(ng))=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic viscosity along geopotential
!  surfaces.  The biharmonic operator is computed by applying
!  the harmonic operator twice.
!-----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients for the BASIC STATE.
!  Notice the recursive storage sequence. It is assumed here that
!  the mixing coefficients are the squared root of the biharmonic
!  viscosity coefficient. For momentum balance purposes, the thickness
!  "Hz" appears only when computing the second harmonic operator.
!  The vertical placement of the gradients is:
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
      K_LOOP1 : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
              cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
              UFx(i,j)=cff*(z_r(i  ,j,k+1)-                             &
     &                      z_r(i-1,j,k+1))
            END DO
          END DO
          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
              cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
              VFe(i,j)=cff*(z_r(i,j  ,k+1)-                             &
     &                      z_r(i,j-1,k+1))
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              dZdx_p(i,j,k2)=0.5_r8*(UFx(i,j-1)+                        &
     &                               UFx(i,j  ))
              dZde_p(i,j,k2)=0.5_r8*(VFe(i-1,j)+                        &
     &                               VFe(i  ,j))
            END DO
          END DO
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              dZdx_r(i,j,k2)=0.5_r8*(UFx(i  ,j)+                        &
     &                               UFx(i+1,j))
              dZde_r(i,j,k2)=0.5_r8*(VFe(i,j  )+                        &
     &                               VFe(i,j+1))
            END DO
          END DO
!
!  Compute BASIC STATE momentum horizontal (1/m/s) and vertical
!  (1/s) gradients.
!
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff=0.5_r8*pm(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
              dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*                 &
     &                           u(i+1,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j)+pn(i  ,j))*                 &
     &                           u(i  ,j,k+1,nrhs))
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                      &
     &                      pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
              dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*               &
     &                           u(i,j  ,k+1,nrhs)-                     &
     &                           (pm(i-1,j-1)+pm(i,j-1))*               &
     &                           u(i,j-1,k+1,nrhs))
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                      &
     &                      pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
              dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*               &
     &                           v(i  ,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j-1)+pn(i-1,j))*               &
     &                           v(i-1,j,k+1,nrhs))
            END DO
          END DO

          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff=0.5_r8*pn(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
              dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*                 &
     &                           v(i,j+1,k+1,nrhs)-                     &
     &                           (pm(i,j-1)+pm(i,j  ))*                 &
     &                           v(i,j  ,k+1,nrhs))
            END DO
          END DO
        END IF

        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
              dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
              dVdz(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstrm1,Jendp1
            DO i=IstrUm1,Iendp1
              UFsx(i,j,k2)=0.0_r8
              UFse(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrVm1,Jendp1
            DO i=Istrm1,Iendp1
              VFsx(i,j,k2)=0.0_r8
              VFse(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-                       &
     &                            z_r(i-1,j,k  )+                       &
     &                            z_r(i  ,j,k+1)-                       &
     &                            z_r(i  ,j,k  )))
              dUdz(i,j,k2)=cff*(u(i,j,k+1,nrhs)-                        &
     &                          u(i,j,k  ,nrhs))
            END DO
          END DO

          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-                       &
     &                            z_r(i,j-1,k  )+                       &
     &                            z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j  ,k  )))
              dVdz(i,j,k2)=cff*(v(i,j,k+1,nrhs)-                        &
     &                          v(i,j,k  ,nrhs))
            END DO
          END DO
        END IF
!
!  Compute BASIC STATE components of the rotated viscous flux
!  (m^4 s-^3/2) along geopotential surfaces in the XI- and
!  ETA-directions.
!
        IF (k.gt.0) THEN
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
              cff3=MIN(dZde_r(i,j,k1),0.0_r8)
              cff4=MAX(dZde_r(i,j,k1),0.0_r8)
              cff=on_r(i,j)*(dnUdx(i,j,k1)-                             &
     &                       0.5_r8*pn(i,j)*                            &
     &                       (cff1*(dUdz(i  ,j,k1)+                     &
     &                              dUdz(i+1,j,k2))+                    &
     &                        cff2*(dUdz(i  ,j,k2)+                     &
     &                              dUdz(i+1,j,k1))))-                  &
     &            om_r(i,j)*(dmVde(i,j,k1)-                             &
     &                       0.5_r8*pm(i,j)*                            &
     &                       (cff3*(dVdz(i,j  ,k1)+                     &
     &                              dVdz(i,j+1,k2))+                    &
     &                        cff4*(dVdz(i,j  ,k2)+                     &
     &                              dVdz(i,j+1,k1))))
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
              UFx(i,j)=on_r(i,j)*on_r(i,j)*Uvis3d_r(i,j,k)*cff
              VFe(i,j)=om_r(i,j)*om_r(i,j)*Vvis3d_r(i,j,k)*cff
# else
              UFx(i,j)=on_r(i,j)*on_r(i,j)*visc3d_r(i,j,k)*cff
              VFe(i,j)=om_r(i,j)*om_r(i,j)*visc3d_r(i,j,k)*cff
# endif
#else
              UFx(i,j)=on_r(i,j)*on_r(i,j)*visc4_r(i,j)*cff
              VFe(i,j)=om_r(i,j)*om_r(i,j)*visc4_r(i,j)*cff
#endif
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              pm_p=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                      &
     &                      pm(i  ,j-1)+pm(i  ,j))
              pn_p=0.25_r8*(pn(i-1,j-1)+pn(i-1,j)+                      &
     &                      pn(i  ,j-1)+pn(i  ,j))
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
              cff3=MIN(dZde_p(i,j,k1),0.0_r8)
              cff4=MAX(dZde_p(i,j,k1),0.0_r8)
              cff=on_p(i,j)*(dnVdx(i,j,k1)-                             &
     &                       0.5_r8*pn_p*                               &
     &                       (cff1*(dVdz(i-1,j,k1)+                     &
     &                              dVdz(i  ,j,k2))+                    &
     &                        cff2*(dVdz(i-1,j,k2)+                     &
     &                              dVdz(i  ,j,k1))))+                  &
     &            om_p(i,j)*(dmUde(i,j,k1)-                             &
     &                       0.5_r8*pm_p*                               &
     &                       (cff3*(dUdz(i,j-1,k1)+                     &
     &                              dUdz(i,j  ,k2))+                    &
     &                        cff4*(dUdz(i,j-1,k2)+                     &
     &                              dUdz(i,j  ,k1))))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
              Uvis_p=0.25_r8*                                           &
     &               (Uvis3d_r(i-1,j-1,k)+Uvis3d_r(i-1,j,k)+            &
     &                Uvis3d_r(i  ,j-1,k)+Uvis3d_r(i  ,j,k))
              Vvis_p=0.25_r8*                                           &
     &               (Vvis3d_r(i-1,j-1,k)+Vvis3d_r(i-1,j,k)+            &
     &                Vvis3d_r(i  ,j-1,k)+Vvis3d_r(i  ,j,k))
              UFe(i,j)=om_p(i,j)*om_p(i,j)*Uvis_p*cff
              VFx(i,j)=on_p(i,j)*on_p(i,j)*Vvis_p*cff
# else
              visc_p=0.25_r8*                                           &
     &               (visc3d_r(i-1,j-1,k)+visc3d_r(i-1,j,k)+            &
     &                visc3d_r(i  ,j-1,k)+visc3d_r(i  ,j,k))
              UFe(i,j)=om_p(i,j)*om_p(i,j)*visc_p*cff
              VFx(i,j)=on_p(i,j)*on_p(i,j)*visc_p*cff
# endif
#else
              UFe(i,j)=om_p(i,j)*om_p(i,j)*visc4_p(i,j)*cff
              VFx(i,j)=on_p(i,j)*on_p(i,j)*visc4_p(i,j)*cff
#endif
            END DO
          END DO
!
!  Compute BASIC STATE vertical flux (m^2 s^-3/2) due to sloping
!  terrain-following surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstrm1,Jendp1
              DO i=IstrUm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Uvis3d_r(i-1,j,k  )+Uvis3d_r(i,j,k  )+             &
     &               Uvis3d_r(i-1,j,k+1)+Uvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i-1,j,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i-1,j,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#else
                cff=0.25_r8*(visc4_r(i-1,j)+visc4_r(i,j))
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#endif
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))

                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                UFsx(i,j,k2)=fac1*                                      &
     &                       (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+        &
     &                        cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+        &
     &                        cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+        &
     &                        cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))

                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                UFse(i,j,k2)=fac2*                                      &
     &                       (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+        &
     &                        cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+        &
     &                        cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+        &
     &                        cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))

                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
                UFsx(i,j,k2)=UFsx(i,j,k2)+                              &
     &                       fac1*                                      &
     &                       (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+        &
     &                        cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+        &
     &                        cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+        &
     &                        cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))

                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
                UFse(i,j,k2)=UFse(i,j,k2)-                              &
     &                       fac2*                                      &
     &                       (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+        &
     &                        cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+        &
     &                        cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+        &
     &                        cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
              END DO
            END DO
!
            DO j=JstrVm1,Jendp1
              DO i=Istrm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Vvis3d_r(i,j-1,k  )+Vvis3d_r(i,j,k  )+             &
     &               Vvis3d_r(i,j-1,k+1)+Vvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i,j-1,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i,j-1,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#else
                cff=0.25_r8*(visc4_r(i,j-1)+visc4_r(i,j))
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#endif
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)

                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                VFsx(i,j,k2)=fac1*                                      &
     &                       (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+        &
     &                        cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+        &
     &                        cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+        &
     &                        cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))

                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                VFse(i,j,k2)=fac2*                                      &
     &                       (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+        &
     &                        cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+        &
     &                        cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+        &
     &                        cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))

                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
                VFsx(i,j,k2)=VFsx(i,j,k2)-                              &
     &                       fac1*                                      &
     &                       (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+        &
     &                        cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+        &
     &                        cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+        &
     &                        cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))

                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
                VFse(i,j,k2)=VFse(i,j,k2)+                              &
     &                       fac2*                                      &
     &                       (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+        &
     &                        cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+        &
     &                        cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+        &
     &                        cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
              END DO
            END DO
          END IF
!
!  Compute BASIC STATE first harmonic operator (m s^-3/2).
!
          DO j=Jstrm1,Jendp1
            DO i=IstrUm1,Iendp1
              cff=0.125_r8*(pm(i-1,j)+pm(i,j))*                         &
     &                     (pn(i-1,j)+pn(i,j))
              cff1=1.0_r8/(0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k)))
              LapU(i,j,k)=cff*((pn(i-1,j)+pn(i,j))*                     &
                               (UFx(i,j)-UFx(i-1,j))+                   &
     &                         (pm(i-1,j)+pm(i,j))*                     &
     &                         (UFe(i,j+1)-UFe(i,j)))+                  &
     &                    cff1*((UFsx(i,j,k2)+UFse(i,j,k2))-            &
     &                          (UFsx(i,j,k1)+UFse(i,j,k1)))
#ifdef MASKING
              LapU(i,j,k)=LapU(i,j,k)*umask(i,j)
#endif
            END DO
          END DO

          DO j=JstrVm1,Jendp1
            DO i=Istrm1,Iendp1
              cff=0.125_r8*(pm(i,j)+pm(i,j-1))*                         &
     &                     (pn(i,j)+pn(i,j-1))
              cff1=1.0_r8/(0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k)))
              LapV(i,j,k)=cff*((pn(i,j-1)+pn(i,j))*                     &
     &                         (VFx(i+1,j)-VFx(i,j))-                   &
     &                         (pm(i,j-1)+pm(i,j))*                     &
     &                         (VFe(i,j)-VFe(i,j-1)))+                  &
     &                    cff1*((VFsx(i,j,k2)+VFse(i,j,k2))-            &
     &                          (VFsx(i,j,k1)+VFse(i,j,k1)))
#ifdef MASKING
              LapV(i,j,k)=LapV(i,j,k)*vmask(i,j)
#endif
            END DO
          END DO
        END IF
      END DO K_LOOP1
!
!  Apply boundary conditions (closed or gradient; except periodic)
!  to the first harmonic operator.
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (ad_LBC(iwest,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                LapU(IstrU-1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                LapU(IstrU-1,j,k)=LapU(IstrU,j,k)
              END DO
            END DO
          END IF
          IF (ad_LBC(iwest,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
                LapV(Istr-1,j,k)=gamma2(ng)*LapV(Istr,j,k)
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
                LapV(Istr-1,j,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (ad_LBC(ieast,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                LapU(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                LapU(Iend+1,j,k)=LapU(Iend,j,k)
              END DO
            END DO
          END IF
          IF (ad_LBC(ieast,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
                LapV(Iend+1,j,k)=gamma2(ng)*LapV(Iend,j,k)
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
                LapV(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (ad_LBC(isouth,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
                LapU(i,Jstr-1,k)=gamma2(ng)*LapU(i,Jstr,k)
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
                LapU(i,Jstr-1,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(isouth,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                LapV(i,JstrV-1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                LapV(i,JstrV-1,k)=LapV(i,JstrV,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (ad_LBC(inorth,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
                LapU(i,Jend+1,k)=gamma2(ng)*LapU(i,Jend,k)
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
                LapU(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(inorth,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                LapV(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                LapV(i,Jend+1,k)=LapV(i,Jend,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            LapU(Istr  ,Jstr-1,k)=0.5_r8*                               &
     &                            (LapU(Istr+1,Jstr-1,k)+               &
     &                             LapU(Istr  ,Jstr  ,k))
            LapV(Istr-1,Jstr  ,k)=0.5_r8*                               &
     &                            (LapV(Istr-1,Jstr+1,k)+               &
     &                             LapV(Istr  ,Jstr  ,k))
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            LapU(Iend+1,Jstr-1,k)=0.5_r8*                               &
     &                            (LapU(Iend  ,Jstr-1,k)+               &
     &                             LapU(Iend+1,Jstr  ,k))
            LapV(Iend+1,Jstr  ,k)=0.5_r8*                               &
     &                            (LapV(Iend  ,Jstr  ,k)+               &
     &                             LapV(Iend+1,Jstr+1,k))
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            LapU(Istr  ,Jend+1,k)=0.5_r8*                               &
     &                            (LapU(Istr+1,Jend+1,k)+               &
     &                             LapU(Istr  ,Jend  ,k))
            LapV(Istr-1,Jend+1,k)=0.5_r8*                               &
     &                            (LapV(Istr  ,Jend+1,k)+               &
     &                             LapV(Istr-1,Jend  ,k))
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            LapU(Iend+1,Jend+1,k)=0.5_r8*                               &
     &                            (LapU(Iend  ,Jend+1,k)+               &
     &                             LapU(Iend+1,Jend  ,k))
            LapV(Iend+1,Jend+1,k)=0.5_r8*                               &
     &                            (LapV(Iend  ,Jend+1,k)+               &
     &                             LapV(Iend+1,Jend  ,k))
          END DO
        END IF
      END IF
!
! Compute adjoint of starting storage recursive indices k1 and k2.
!
      k1=2
      k2=1
      DO k=0,N(ng)
        k1=k2
        k2=3-k1
      END DO
!
!  Compute required BASIC STATE fields. Need to look forward in
!  recursive kk index.
!
      K_LOOP2 : DO k=N(ng),0,-1
        k2b=1
        DO kk=0,k
          k1b=k2b
          k2b=3-k1b
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          IF (kk.lt.N(ng)) THEN
            DO j=Jstrm2,Jendp2
              DO i=IstrUm2,Iendp2
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
                UFx(i,j)=cff*(z_r(i  ,j,kk+1)-                          &
     &                        z_r(i-1,j,kk+1))
              END DO
            END DO
            DO j=JstrVm2,Jendp2
              DO i=Istrm2,Iendp2
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
                VFe(i,j)=cff*(z_r(i,j  ,kk+1)-                          &
     &                        z_r(i,j-1,kk+1))
              END DO
            END DO
!
            DO j=Jstrm1,Jendp2
               DO i=Istrm1,Iendp2
                dZdx_p(i,j,k2b)=0.5_r8*(UFx(i,j-1)+                     &
     &                                  UFx(i,j  ))
                dZde_p(i,j,k2b)=0.5_r8*(VFe(i-1,j)+                     &
     &                                  VFe(i  ,j))
              END DO
            END DO
            DO j=JstrVm2,Jendp1
              DO i=IstrUm2,Iendp1
                dZdx_r(i,j,k2b)=0.5_r8*(UFx(i  ,j)+                     &
     &                                  UFx(i+1,j))
                dZde_r(i,j,k2b)=0.5_r8*(VFe(i,j  )+                     &
     &                                  VFe(i,j+1))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=Jstrm1,Jendp2
                DO i=Istrm1,Iendp2
                  dZdx_p(i,j,k1b)=0.0_r8
                  dZde_p(i,j,k1b)=0.0_r8
                END DO
              END DO
              DO j=JstrVm2,Jendp1
                DO i=IstrUm2,Iendp1
                  dZdx_r(i,j,k1b)=0.0_r8
                  dZde_r(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
!
!  Compute momentum horizontal (m^-1 s^-3/2) and vertical (s^-3/2)
!  gradients.
!
            DO j=JstrV-1,Jend
              DO i=IstrU-1,Iend
                cff=0.5_r8*pm(i,j)
#ifdef MASKING
                cff=cff*rmask(i,j)
#endif
                dnUdx(i,j,k2b)=cff*((pn(i  ,j)+pn(i+1,j))*              &
     &                              LapU(i+1,j,kk+1)-                   &
     &                              (pn(i-1,j)+pn(i  ,j))*              &
     &                              LapU(i  ,j,kk+1))
              END DO
            END DO

            DO j=Jstr,Jend+1
              DO i=Istr,Iend+1
                cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                    &
     &                        pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
                cff=cff*pmask(i,j)
#endif
                dmUde(i,j,k2b)=cff*((pm(i-1,j  )+pm(i,j  ))*            &
     &                              LapU(i,j  ,kk+1)-                   &
     &                              (pm(i-1,j-1)+pm(i,j-1))*            &
     &                              LapU(i,j-1,kk+1))
              END DO
            END DO

            DO j=Jstr,Jend+1
              DO i=Istr,Iend+1
                cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                    &
     &                        pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
                cff=cff*pmask(i,j)
#endif
                dnVdx(i,j,k2b)=cff*((pn(i  ,j-1)+pn(i  ,j))*            &
     &                              LapV(i  ,j,kk+1)-                   &
     &                              (pn(i-1,j-1)+pn(i-1,j))*            &
     &                              LapV(i-1,j,kk+1))
              END DO
            END DO

            DO j=JstrV-1,Jend
              DO i=IstrU-1,Iend
                cff=0.5_r8*pn(i,j)
#ifdef MASKING
                cff=cff*rmask(i,j)
#endif
                dmVde(i,j,k2b)=cff*((pm(i,j  )+pm(i,j+1))*              &
     &                              LapV(i,j+1,kk+1)-                   &
     &                              (pm(i,j-1)+pm(i,j  ))*              &
     &                              LapV(i,j  ,kk+1))
              END DO
            END DO

            IF (kk.eq.0) THEN
              DO j=JstrV-1,Jend
                DO i=IstrU-1,Iend
                  dnUdx(i,j,k1b)=0.0_r8
                END DO
              END DO
              DO j=Jstr,Jend+1
                DO i=Istr,Iend+1
                  dmUde(i,j,k1b)=0.0_r8
                END DO
              END DO
              DO j=Jstr,Jend+1
                DO i=Istr,Iend+1
                  dnVdx(i,j,k1b)=0.0_r8
                END DO
              END DO
              DO j=JstrV-1,Jend
                DO i=IstrU-1,Iend
                  dmVde(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
          END IF

          IF ((kk.eq.0).or.(kk.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                dUdz(i,j,k2b)=0.0_r8
              END DO
            END DO
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                dVdz(i,j,k2b)=0.0_r8
              END DO
            END DO

            IF (kk.eq.0) THEN
              DO j=Jstr-1,Jend+1
                DO i=IstrU-1,Iend+1
                  dUdz(i,j,k1b)=0.0_r8
                END DO
              END DO
              DO j=JstrV-1,Jend+1
                DO i=Istr-1,Iend+1
                  dVdz(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                cff=1.0_r8/(0.5_r8*(z_r(i-1,j,kk+1)-                    &
     &                              z_r(i-1,j,kk  )+                    &
     &                              z_r(i  ,j,kk+1)-                    &
     &                              z_r(i  ,j,kk  )))
                dUdz(i,j,k2b)=cff*(LapU(i,j,kk+1)-                      &
     &                             LapU(i,j,kk  ))
              END DO
            END DO

            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(0.5_r8*(z_r(i,j-1,kk+1)-                    &
     &                              z_r(i,j-1,kk  )+                    &
     &                              z_r(i,j  ,kk+1)-                    &
     &                              z_r(i,j  ,kk  )))
                dVdz(i,j,k2b)=cff*(LapV(i,j,kk+1)-                      &
     &                             LapV(i,j,kk  ))
              END DO
            END DO
          END IF
        END DO
!
        IF (k.gt.0) THEN
!
!  Time-step biharmonic, geopotential viscosity term. Notice that
!  momentum at this stage is HzU and HzV and has m2/s units.  Add
!  contribution for barotropic forcing terms.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=dt(ng)*0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!#ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3yvis)= cff*cff2-dt(ng)*cff4
!!            DiaV3wrk(i,j,k,M3xvis)=-cff*cff1-dt(ng)*cff3
!!            DiaV3wrk(i,j,k,M3hvis)=-cff5-cff6
!!            DiaRVfrc(i,j,3,M2yvis)=DiaRVfrc(i,j,3,M2yvis)+cff2-cff4
!!            DiaRVfrc(i,j,3,M2xvis)=DiaRVfrc(i,j,3,M2xvis)-cff1-cff3
!!            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-cff1+cff2-  &
!!   &                                                      cff3-cff4
!!#endif
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_cff5-tl_cff6
!>
              ad_cff5=ad_cff5-ad_v(i,j,k,nnew)
              ad_cff6=ad_cff6-ad_v(i,j,k,nnew)
!>            tl_rvfrc(i,j)=tl_rvfrc(i,j)-                              &
!>   &                      tl_cff1+tl_cff2-tl_cff3-tl_cff4
!>
              ad_cff1=ad_cff1-ad_rvfrc(i,j)
              ad_cff2=ad_cff2+ad_rvfrc(i,j)
              ad_cff3=ad_cff3-ad_rvfrc(i,j)
              ad_cff4=ad_cff4-ad_rvfrc(i,j)
!>            tl_cff6=dt(ng)*(tl_cff3+tl_cff4)
!>
              adfac=dt(ng)*ad_cff6
              ad_cff3=ad_cff3+adfac
              ad_cff4=ad_cff4+adfac
              ad_cff6=0.0_r8
!>            tl_cff5=cff*(tl_cff1-tl_cff2)
!>
              adfac=cff*ad_cff5
              ad_cff1=ad_cff1+adfac
              ad_cff2=ad_cff2-adfac
              ad_cff5=0.0_r8
!>            tl_cff4=tl_VFse(i,j,k2)-tl_VFse(i,j,k1)
!>
              ad_VFse(i,j,k1)=ad_VFse(i,j,k1)-ad_cff4
              ad_VFse(i,j,k2)=ad_VFse(i,j,k2)+ad_cff4
              ad_cff4=0.0_r8
!>            tl_cff3=tl_VFsx(i,j,k2)-tl_VFsx(i,j,k1)
!>
              ad_VFsx(i,j,k1)=ad_VFsx(i,j,k1)-ad_cff3
              ad_VFsx(i,j,k2)=ad_VFsx(i,j,k2)+ad_cff3
              ad_cff3=0.0_r8
!>            tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                       &
!>   &                (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>
              adfac=0.5_r8*(pm(i,j-1)+pm(i,j))*ad_cff2
              ad_VFe(i,j-1)=ad_VFe(i,j-1)-adfac
              ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
              ad_cff2=0.0_r8
!>            tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                       &
!>   &                (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>
              adfac=0.5_r8*(pn(i,j-1)+pn(i,j))*ad_cff1
              ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac
              ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac
              ad_cff1=0.0_r8
            END DO
          END DO

          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=dt(ng)*0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
#ifdef DIAGNOSTICS_UV
!!            DiaU3wrk(i,j,k,M3yvis)=-cff*cff2-dt(ng)*cff4
!!            DiaU3wrk(i,j,k,M3xvis)=-cff*cff1-dt(ng)*cff3
!!            DiaU3wrk(i,j,k,M3hvis)=-cff5-cff6
!!            DiaRUfrc(i,j,3,M2yvis)=DiaRUfrc(i,j,3,M2yvis)-cff2-cff4
!!            DiaRUfrc(i,j,3,M2xvis)=DiaRUfrc(i,j,3,M2xvis)-cff1-cff3
!!            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-cff1-cff2-  &
!!   &                                                      cff3-cff4
#endif
!>            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_cff5-tl_cff6
!>
              ad_cff5=ad_cff5-ad_u(i,j,k,nnew)
              ad_cff6=ad_cff6-ad_u(i,j,k,nnew)
!>            tl_rufrc(i,j)=tl_rufrc(i,j)-                              &
!>   &                      tl_cff1-tl_cff2-tl_cff3-tl_cff4
!>
              ad_cff1=ad_cff1-ad_rufrc(i,j)
              ad_cff2=ad_cff2-ad_rufrc(i,j)
              ad_cff3=ad_cff3-ad_rufrc(i,j)
              ad_cff4=ad_cff4-ad_rufrc(i,j)
!>            tl_cff6=dt(ng)*(tl_cff3+tl_cff4)
!>
              adfac=dt(ng)*ad_cff6
              ad_cff3=ad_cff3+adfac
              ad_cff4=ad_cff4+adfac
              ad_cff6=0.0_r8
!>            tl_cff5=cff*(tl_cff1+tl_cff2)
!>
              adfac=cff*ad_cff5
              ad_cff1=ad_cff1+adfac
              ad_cff2=ad_cff2+adfac
              ad_cff5=0.0_r8
!>            tl_cff4=tl_UFse(i,j,k2)-tl_UFse(i,j,k1)
!>
              ad_UFse(i,j,k1)=ad_UFse(i,j,k1)-ad_cff4
              ad_UFse(i,j,k2)=ad_UFse(i,j,k2)+ad_cff4
              ad_cff4=0.0_r8
!>            tl_cff3=tl_UFsx(i,j,k2)-tl_UFsx(i,j,k1)
!>
              ad_UFsx(i,j,k1)=ad_UFsx(i,j,k1)-ad_cff3
              ad_UFsx(i,j,k2)=ad_UFsx(i,j,k2)+ad_cff3
              ad_cff3=0.0_r8
!>            tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                       &
!>   &                (tl_UFe(i,j+1)-tl_UFe(i  ,j))
!>
              adfac=0.5_r8*(pm(i-1,j)+pm(i,j))*ad_cff2
              ad_UFe(i,j  )=ad_UFe(i,j  )-adfac
              ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac
              ad_cff2=0.0_r8
!>            tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                       &
!>   &                (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>
              adfac=0.5_r8*(pn(i-1,j)+pn(i,j))*ad_cff1
              ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac
              ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
              ad_cff1=0.0_r8
            END DO
          END DO
!
!  Compute vertical flux (m^2 s^-3/2) due to sloping terrain-following
!  surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Vvis3d_r(i,j-1,k  )+Vvis3d_r(i,j,k  )+             &
     &               Vvis3d_r(i,j-1,k+1)+Vvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i,j-1,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i,j-1,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#else
                cff=0.25_r8*(visc4_r(i,j-1)+visc4_r(i,j))
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#endif
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)
!
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+     &
!>   &                           cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+     &
!>   &                           cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+     &
!>   &                           cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+             &
     &                   cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+             &
     &                   cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+             &
     &                   cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))*            &
     &                  ad_VFse(i,j,k2)
#endif
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          fac2*                                   &
!>   &                          (tl_cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+  &
!>   &                           tl_cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+  &
!>   &                           tl_cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+  &
!>   &                           tl_cff4*(cff8*dmUdz-dmUde(i+1,j,k1))+  &
!>   &                           cff1*(tl_cff5*dmUdz+cff5*tl_dmUdz-     &
!>   &                                 tl_dmUde(i  ,j,k1))+             &
!>   &                           cff2*(tl_cff6*dmUdz+cff6*tl_dmUdz-     &
!>   &                                 tl_dmUde(i+1,j,k2))+             &
!>   &                           cff3*(tl_cff7*dmUdz+cff7*tl_dmUdz-     &
!>   &                                 tl_dmUde(i  ,j,k2))+             &
!>   &                           cff4*(tl_cff8*dmUdz+cff8*tl_dmUdz-     &
!>   &                                 tl_dmUde(i+1,j,k1)))
!>
                adfac=fac2*ad_VFse(i,j,k2)
                adfac1=adfac*dmUdz
                ad_cff1=ad_cff1+(cff5*dmUdz-dmUde(i  ,j,k1))*adfac
                ad_cff2=ad_cff2+(cff6*dmUdz-dmUde(i+1,j,k2))*adfac
                ad_cff3=ad_cff3+(cff7*dmUdz-dmUde(i  ,j,k2))*adfac
                ad_cff4=ad_cff4+(cff8*dmUdz-dmUde(i+1,j,k1))*adfac
                ad_cff5=ad_cff5+cff1*adfac1
                ad_cff6=ad_cff6+cff2*adfac1
                ad_cff7=ad_cff7+cff3*adfac1
                ad_cff8=ad_cff8+cff4*adfac1
                ad_dmUdz=ad_dmUdz+                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dmUde(i  ,j,k1)=ad_dmUde(i  ,j,k1)-cff1*adfac
                ad_dmUde(i+1,j,k2)=ad_dmUde(i+1,j,k2)-cff2*adfac
                ad_dmUde(i  ,j,k2)=ad_dmUde(i  ,j,k2)-cff3*adfac
                ad_dmUde(i+1,j,k1)=ad_dmUde(i+1,j,k1)-cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_p(i+1,j,k1)))*        &
!>   &                  tl_dZde_p(i+1,j,k1)
!>
                ad_dZde_p(i+1,j,k1)=ad_dZde_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i+1,j,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_p(i  ,j,k2)))*        &
!>   &                  tl_dZde_p(i  ,j,k2)
!>
                ad_dZde_p(i  ,j,k2)=ad_dZde_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i  ,j,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_p(i+1,j,k2)))*        &
!>   &                  tl_dZde_p(i+1,j,k2)
!>
                ad_dZde_p(i+1,j,k2)=ad_dZde_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i+1,j,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_p(i  ,j,k1)))*        &
!>   &                  tl_dZde_p(i  ,j,k1)
!>
                ad_dZde_p(i  ,j,k1)=ad_dZde_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i  ,j,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
!>   &                  tl_dZdx_p(i+1,j,k1)
!>
                ad_dZdx_p(i+1,j,k1)=ad_dZdx_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
!>   &                  tl_dZdx_p(i  ,j,k2)
!>
                ad_dZdx_p(i  ,j,k2)=ad_dZdx_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
!>   &                  tl_dZdx_p(i+1,j,k2)
!>
                ad_dZdx_p(i+1,j,k2)=ad_dZdx_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
!>   &                  tl_dZdx_p(i  ,j,k1)
!>
                ad_dZdx_p(i  ,j,k1)=ad_dZdx_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8

                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+     &
!>   &                           cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+     &
!>   &                           cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+     &
!>   &                           cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
!>
                ad_fac1=ad_fac1-                                        &
     &                  (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+             &
     &                   cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+             &
     &                   cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+             &
     &                   cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))*            &
     &                  ad_VFsx(i,j,k2)
#endif
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
!>   &                          fac1*                                   &
!>   &                          (tl_cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+  &
!>   &                           tl_cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+  &
!>   &                           tl_cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+  &
!>   &                           tl_cff4*(cff8*dnUdz-dnUdx(i,j  ,k1))+  &
!>   &                           cff1*(tl_cff5*dnUdz+cff5*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j-1,k1))+             &
!>   &                           cff2*(tl_cff6*dnUdz+cff6*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j  ,k2))+             &
!>   &                           cff3*(tl_cff7*dnUdz+cff7*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j-1,k2))+             &
!>   &                           cff4*(tl_cff8*dnUdz+cff8*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j  ,k1)))
!>
                adfac=fac1*ad_VFsx(i,j,k2)
                adfac1=adfac*dnUdz
                ad_cff1=ad_cff1-(cff5*dnUdz-dnUdx(i,j-1,k1))*adfac
                ad_cff2=ad_cff2-(cff6*dnUdz-dnUdx(i,j  ,k2))*adfac
                ad_cff3=ad_cff3-(cff7*dnUdz-dnUdx(i,j-1,k2))*adfac
                ad_cff4=ad_cff4-(cff8*dnUdz-dnUdx(i,j  ,k1))*adfac
                ad_cff5=ad_cff5-cff1*adfac1
                ad_cff6=ad_cff6-cff2*adfac1
                ad_cff7=ad_cff7-cff3*adfac1
                ad_cff8=ad_cff8-cff4*adfac1
                ad_dnUdz=ad_dnUdz-                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dnUdx(i,j-1,k1)=ad_dnUdx(i,j-1,k1)+cff1*adfac
                ad_dnUdx(i,j  ,k2)=ad_dnUdx(i,j  ,k2)+cff2*adfac
                ad_dnUdx(i,j-1,k2)=ad_dnUdx(i,j-1,k2)+cff3*adfac
                ad_dnUdx(i,j  ,k1)=ad_dnUdx(i,j  ,k1)+cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j  ,k1)))*        &
!>   &                  tl_dZdx_r(i,j  ,k1)
!>
                ad_dZdx_r(i,j  ,k1)=ad_dZdx_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i,j  ,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j-1,k2)))*        &
!>   &                  tl_dZdx_r(i,j-1,k2)
!>
                ad_dZdx_r(i,j-1,k2)=ad_dZdx_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i,j-1,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j  ,k2)))*        &
!>   &                  tl_dZdx_r(i,j  ,k2)
!>
                ad_dZdx_r(i,j  ,k2)=ad_dZdx_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i,j  ,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j-1,k1)))*        &
!>   &                  tl_dZdx_r(i,j-1,k1)
!>
                ad_dZdx_r(i,j-1,k1)=ad_dZdx_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i,j-1,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
!>   &                  tl_dZde_r(i,j  ,k1)
!>
                ad_dZde_r(i,j  ,k1)=ad_dZde_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j  ,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
!>   &                  tl_dZde_r(i,j-1,k2)
!>
                ad_dZde_r(i,j-1,k2)=ad_dZde_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j-1,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
!>   &                  tl_dZde_r(i,j  ,k2)
!>
                ad_dZde_r(i,j  ,k2)=ad_dZde_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
!>   &                  tl_dZde_r(i,j-1,k1)
!>
                ad_dZde_r(i,j-1,k1)=ad_dZde_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+     &
!>   &                           cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+     &
!>   &                           cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+     &
!>   &                           cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+             &
     &                   cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+             &
     &                   cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+             &
     &                   cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))*            &
     &                  ad_VFse(i,j,k2)
#endif
!>              tl_VFse(i,j,k2)=fac2*                                   &
!>   &                          (tl_cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+  &
!>   &                           tl_cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+  &
!>   &                           tl_cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+  &
!>   &                           tl_cff4*(cff4*dmVdz-dmVde(i,j  ,k1))+  &
!>   &                           cff1*(tl_cff1*dmVdz+cff1*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j-1,k1))+             &
!>   &                           cff2*(tl_cff2*dmVdz+cff2*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j  ,k2))+             &
!>   &                           cff3*(tl_cff3*dmVdz+cff3*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j-1,k2))+             &
!>   &                           cff4*(tl_cff4*dmVdz+cff4*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j  ,k1)))
!>
                cff=2.0_r8*dmVdz
                adfac=fac2*ad_VFse(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dmVde(i,j-1,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dmVde(i,j  ,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dmVde(i,j-1,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dmVde(i,j  ,k1))*adfac
                ad_dmVdz=ad_dmVdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dmVde(i,j-1,k1)=ad_dmVde(i,j-1,k1)-cff1*adfac
                ad_dmVde(i,j  ,k2)=ad_dmVde(i,j  ,k2)-cff2*adfac
                ad_dmVde(i,j-1,k2)=ad_dmVde(i,j-1,k2)-cff3*adfac
                ad_dmVde(i,j  ,k1)=ad_dmVde(i,j  ,k1)-cff4*adfac
                ad_VFse(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
!>   &                  tl_dZde_r(i,j  ,k1)
!>
                ad_dZde_r(i,j  ,k1)=ad_dZde_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j  ,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
!>   &                  tl_dZde_r(i,j-1,k2)
!>
                ad_dZde_r(i,j-1,k2)=ad_dZde_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j-1,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
!>   &                  tl_dZde_r(i,j  ,k2)
!>
                ad_dZde_r(i,j  ,k2)=ad_dZde_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
!>   &                  tl_dZde_r(i,j-1,k1)
!>
                ad_dZde_r(i,j-1,k1)=ad_dZde_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+     &
!>   &                           cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+     &
!>   &                           cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+     &
!>   &                           cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+             &
     &                   cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+             &
     &                   cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+             &
     &                   cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))*            &
     &                  ad_VFsx(i,j,k2)
#endif
!>              tl_VFsx(i,j,k2)=fac1*                                   &
!>   &                          (tl_cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+  &
!>   &                           tl_cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+  &
!>   &                           tl_cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+  &
!>   &                           tl_cff4*(cff4*dnVdz-dnVdx(i+1,j,k1))+  &
!>   &                           cff1*(tl_cff1*dnVdz+cff1*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i  ,j,k1))+             &
!>   &                           cff2*(tl_cff2*dnVdz+cff2*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i+1,j,k2))+             &
!>   &                           cff3*(tl_cff3*dnVdz+cff3*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i  ,j,k2))+             &
!>   &                           cff4*(tl_cff4*dnVdz+cff4*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i+1,j,k1)))
!>
                cff=2.0_r8*dnVdz
                adfac=fac1*ad_VFsx(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dnVdx(i  ,j,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dnVdx(i+1,j,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dnVdx(i  ,j,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dnVdx(i+1,j,k1))*adfac
                ad_dnVdz=ad_dnVdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dnVdx(i  ,j,k1)=ad_dnVdx(i  ,j,k1)-cff1*adfac
                ad_dnVdx(i+1,j,k2)=ad_dnVdx(i+1,j,k2)-cff2*adfac
                ad_dnVdx(i  ,j,k2)=ad_dnVdx(i  ,j,k2)-cff3*adfac
                ad_dnVdx(i+1,j,k1)=ad_dnVdx(i+1,j,k1)-cff4*adfac
                ad_VFsx(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
!>   &                  tl_dZdx_p(i+1,j,k1)
!>
                ad_dZdx_p(i+1,j,k1)=ad_dZdx_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
!>   &                  tl_dZdx_p(i  ,j,k2)
!>
                ad_dZdx_p(i  ,j,k2)=ad_dZdx_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
!>   &                  tl_dZdx_p(i+1,j,k2)
!>
                ad_dZdx_p(i+1,j,k2)=ad_dZdx_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
!>   &                  tl_dZdx_p(i  ,j,k1)
!>
                ad_dZdx_p(i  ,j,k1)=ad_dZdx_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
!>              tl_dmVdz=cff*tl_dVdz(i,j,k2)
!>
                ad_dVdz(i,j,k2)=ad_dVdz(i,j,k2)+cff*ad_dmVdz
                ad_dmVdz=0.0_r8
!>              tl_dmUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
!>   &                                tl_dUdz(i+1,j  ,k2)+              &
!>   &                                tl_dUdz(i  ,j-1,k2)+              &
!>   &                                tl_dUdz(i+1,j-1,k2))
!>
                adfac=cff*0.25_r8*ad_dmUdz
                ad_dUdz(i  ,j-1,k2)=ad_dUdz(i  ,j-1,k2)+adfac
                ad_dUdz(i+1,j-1,k2)=ad_dUdz(i+1,j-1,k2)+adfac
                ad_dUdz(i  ,j  ,k2)=ad_dUdz(i  ,j  ,k2)+adfac
                ad_dUdz(i+1,j  ,k2)=ad_dUdz(i+1,j  ,k2)+adfac
                ad_dmUdz=0.0_r8
!
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
!>              tl_dnVdz=cff*tl_dVdz(i,j,k2)
!>
                ad_dVdz(i,j,k2)=ad_dVdz(i,j,k2)+cff*ad_dnVdz
                ad_dnVdz=0.0_r8
!>              tl_dnUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
!>   &                                tl_dUdz(i+1,j  ,k2)+              &
!>   &                                tl_dUdz(i  ,j-1,k2)+              &
!>   &                                tl_dUdz(i+1,j-1,k2))
!>
                adfac=cff*0.25_r8*ad_dnUdz
                ad_dUdz(i  ,j-1,k2)=ad_dUdz(i  ,j-1,k2)+adfac
                ad_dUdz(i+1,j-1,k2)=ad_dUdz(i+1,j-1,k2)+adfac
                ad_dUdz(i  ,j  ,k2)=ad_dUdz(i  ,j  ,k2)+adfac
                ad_dUdz(i+1,j  ,k2)=ad_dUdz(i+1,j  ,k2)+adfac
                ad_dnUdz=0.0_r8
#ifdef VISC_3DCOEF
!>              tl_fac2=tl_cff*om_v(i,j)
!>              tl_fac1=tl_cff*on_v(i,j)
!>
                ad_cff=ad_cff+                                          &
     &                 on_v(i,j)*ad_fac1+om_v(i,j)*ad_fac2
                ad_fac1=0.0_r8
                ad_fac2=0.0_r8
# ifdef UV_U3ADV_SPLIT
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_Vvis3d_r(i,j-1,k  )+tl_Vvis3d_r(i,j,k  )+    &
!>   &                  tl_Vvis3d_r(i,j-1,k+1)+tl_Vvis3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_Vvis3d_r(i,j-1,k  )=ad_Vvis3d_r(i,j-1,k  )+adfac
                ad_Vvis3d_r(i,j  ,k  )=ad_Vvis3d_r(i,j  ,k  )+adfac
                ad_Vvis3d_r(i,j-1,k+1)=ad_Vvis3d_r(i,j-1,k+1)+adfac
                ad_Vvis3d_r(i,j  ,k+1)=ad_Vvis3d_r(i,j  ,k+1)+adfac
                ad_cff=0.0_r8
# else
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_visc3d_r(i,j-1,k  )+tl_visc3d_r(i,j,k  )+    &
!>   &                  tl_visc3d_r(i,j-1,k+1)+tl_visc3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_visc3d_r(i,j-1,k  )=ad_visc3d_r(i,j-1,k  )+adfac
                ad_visc3d_r(i,j  ,k  )=ad_visc3d_r(i,j  ,k  )+adfac
                ad_visc3d_r(i,j-1,k+1)=ad_visc3d_r(i,j-1,k+1)+adfac
                ad_visc3d_r(i,j  ,k+1)=ad_visc3d_r(i,j  ,k+1)+adfac
                ad_cff=0.0_r8
# endif
#endif
              END DO
            END DO
!
            DO j=Jstr,Jend
              DO i=IstrU,Iend
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Uvis3d_r(i-1,j,k  )+Uvis3d_r(i,j,k  )+             &
     &               Uvis3d_r(i-1,j,k+1)+Uvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i-1,j,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i-1,j,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#else
                cff=0.25_r8*(visc4_r(i-1,j)+visc4_r(i,j))
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#endif
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))

                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+     &
!>   &                           cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+     &
!>   &                           cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+     &
!>   &                           cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
!>
                ad_fac2=ad_fac2-                                        &
     &                  (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+             &
     &                   cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+             &
     &                   cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+             &
     &                   cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))*            &
     &                  ad_UFse(i,j,k2)
#endif
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
!>   &                          fac2*                                   &
!>   &                          (tl_cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+  &
!>   &                           tl_cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+  &
!>   &                           tl_cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+  &
!>   &                           tl_cff4*(cff8*dmVdz-dmVde(i  ,j,k1))+  &
!>   &                           cff1*(tl_cff5*dmVdz+cff5*tl_dmVdz-     &
!>   &                                 tl_dmVde(i-1,j,k1))+             &
!>   &                           cff2*(tl_cff6*dmVdz+cff6*tl_dmVdz-     &
!>   &                                 tl_dmVde(i  ,j,k2))+             &
!>   &                           cff3*(tl_cff7*dmVdz+cff7*tl_dmVdz-     &
!>   &                                 tl_dmVde(i-1,j,k2))+             &
!>   &                           cff4*(tl_cff8*dmVdz+cff8*tl_dmVdz-     &
!>   &                                 tl_dmVde(i  ,j,k1)))
!>
                adfac=fac2*ad_UFse(i,j,k2)
                adfac1=adfac*dmVdz
                ad_cff1=ad_cff1-(cff5*dmVdz-dmVde(i-1,j,k1))*adfac
                ad_cff2=ad_cff2-(cff6*dmVdz-dmVde(i  ,j,k2))*adfac
                ad_cff3=ad_cff3-(cff7*dmVdz-dmVde(i-1,j,k2))*adfac
                ad_cff4=ad_cff4-(cff8*dmVdz-dmVde(i  ,j,k1))*adfac
                ad_cff5=ad_cff5-cff1*adfac1
                ad_cff6=ad_cff6-cff2*adfac1
                ad_cff7=ad_cff7-cff3*adfac1
                ad_cff8=ad_cff8-cff4*adfac1
                ad_dmVdz=ad_dmVdz-                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dmVde(i-1,j,k1)=ad_dmVde(i-1,j,k1)+cff1*adfac
                ad_dmVde(i  ,j,k2)=ad_dmVde(i  ,j,k2)+cff2*adfac
                ad_dmVde(i-1,j,k2)=ad_dmVde(i-1,j,k2)+cff3*adfac
                ad_dmVde(i  ,j,k1)=ad_dmVde(i  ,j,k1)+cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_r(i  ,j,k1)))*        &
!>   &                  tl_dZde_r(i  ,j,k1)
!>
                ad_dZde_r(i  ,j,k1)=ad_dZde_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i  ,j,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_r(i-1,j,k2)))*        &
!>   &                  tl_dZde_r(i-1,j,k2)
!>
                ad_dZde_r(i-1,j,k2)=ad_dZde_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i-1,j,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_r(i  ,j,k2)))*        &
!>   &                  tl_dZde_r(i  ,j,k2)
!>
                ad_dZde_r(i  ,j,k2)=ad_dZde_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i  ,j,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_r(i-1,j,k1)))*        &
!>   &                  tl_dZde_r(i-1,j,k1)
!>
                ad_dZde_r(i-1,j,k1)=ad_dZde_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i-1,j,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
!>   &                  tl_dZdx_r(i  ,j,k1)
!>
                ad_dZdx_r(i  ,j,k1)=ad_dZdx_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
!>   &                  tl_dZdx_r(i-1,j,k2)
!>
                ad_dZdx_r(i-1,j,k2)=ad_dZdx_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
!>   &                  tl_dZdx_r(i  ,j,k2)
!>
                ad_dZdx_r(i  ,j,k2)=ad_dZdx_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
!>   &                  tl_dZdx_r(i-1,j,k1)
!>
                ad_dZdx_r(i-1,j,k1)=ad_dZdx_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+     &
!>   &                           cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+     &
!>   &                           cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+     &
!>   &                           cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+             &
     &                   cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+             &
     &                   cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+             &
     &                   cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))*            &
     &                  ad_UFsx(i,j,k2)
#endif
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          fac1*                                   &
!>   &                          (tl_cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+  &
!>   &                           tl_cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+  &
!>   &                           tl_cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+  &
!>   &                           tl_cff4*(cff8*dnVdz-dnVdx(i,j+1,k1))+  &
!>   &                           cff1*(tl_cff5*dnVdz+cff5*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j  ,k1))+             &
!>   &                           cff2*(tl_cff6*dnVdz+cff6*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j+1,k2))+             &
!>   &                           cff3*(tl_cff7*dnVdz+cff7*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j  ,k2))+             &
!>   &                           cff4*(tl_cff8*dnVdz+cff8*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j+1,k1)))
!>
                adfac=fac1*ad_UFsx(i,j,k2)
                adfac1=adfac*dnVdz
                ad_cff1=ad_cff1+(cff5*dnVdz-dnVdx(i,j  ,k1))*adfac
                ad_cff2=ad_cff2+(cff6*dnVdz-dnVdx(i,j+1,k2))*adfac
                ad_cff3=ad_cff3+(cff7*dnVdz-dnVdx(i,j  ,k2))*adfac
                ad_cff4=ad_cff4+(cff8*dnVdz-dnVdx(i,j+1,k1))*adfac
                ad_cff5=ad_cff5+cff1*adfac1
                ad_cff6=ad_cff6+cff2*adfac1
                ad_cff7=ad_cff7+cff3*adfac1
                ad_cff8=ad_cff8+cff4*adfac1
                ad_dnVdz=ad_dnVdz+                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dnVdx(i,j  ,k1)=ad_dnVdx(i,j  ,k1)-cff1*adfac
                ad_dnVdx(i,j+1,k2)=ad_dnVdx(i,j+1,k2)-cff2*adfac
                ad_dnVdx(i,j  ,k2)=ad_dnVdx(i,j  ,k2)-cff3*adfac
                ad_dnVdx(i,j+1,k1)=ad_dnVdx(i,j+1,k1)-cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j+1,k1)))*        &
!>   &                  tl_dZdx_p(i,j+1,k1)
!>
                ad_dZdx_p(i,j+1,k1)=ad_dZdx_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i,j+1,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j  ,k2)))*        &
!>   &                  tl_dZdx_p(i,j  ,k2)
!>
                ad_dZdx_p(i,j  ,k2)=ad_dZdx_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i,j  ,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j+1,k2)))*        &
!>   &                  tl_dZdx_p(i,j+1,k2)
!>
                ad_dZdx_p(i,j+1,k2)=ad_dZdx_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i,j+1,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j  ,k1)))*        &
!>   &                  tl_dZdx_p(i,j  ,k1)
!>
                ad_dZdx_p(i,j  ,k1)=ad_dZdx_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i,j  ,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
!>   &                  tl_dZde_p(i,j+1,k1)
!>
                ad_dZde_p(i,j+1,k1)=ad_dZde_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j+1,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
!>   &                  tl_dZde_p(i,j  ,k2)
!>
                ad_dZde_p(i,j  ,k2)=ad_dZde_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j  ,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
!>   &                  tl_dZde_p(i,j+1,k2)
!>
                ad_dZde_p(i,j+1,k2)=ad_dZde_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
!>   &                  tl_dZde_p(i,j  ,k1)
!>
                ad_dZde_p(i,j  ,k1)=ad_dZde_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)+
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+     &
!>   &                           cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+     &
!>   &                           cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+     &
!>   &                           cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+             &
     &                   cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+             &
     &                   cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+             &
     &                   cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))*            &
     &                  ad_UFse(i,j,k2)
#endif
!>              tl_UFse(i,j,k2)=fac2*                                   &
!>   &                          (tl_cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+  &
!>   &                           tl_cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+  &
!>   &                           tl_cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+  &
!>   &                           tl_cff4*(cff4*dmUdz-dmUde(i,j+1,k1))+  &
!>   &                           cff1*(tl_cff1*dmUdz+cff1*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j  ,k1))+             &
!>   &                           cff2*(tl_cff2*dmUdz+cff2*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j+1,k2))+             &
!>   &                           cff3*(tl_cff3*dmUdz+cff3*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j  ,k2))+             &
!>   &                           cff4*(tl_cff4*dmUdz+cff4*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j+1,k1)))
!>
                cff=2.0_r8*dmUdz
                adfac=fac2*ad_UFse(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dmUde(i,j  ,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dmUde(i,j+1,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dmUde(i,j  ,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dmUde(i,j+1,k1))*adfac
                ad_dmUdz=ad_dmUdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dmUde(i,j  ,k1)=ad_dmUde(i,j  ,k1)-cff1*adfac
                ad_dmUde(i,j+1,k2)=ad_dmUde(i,j+1,k2)-cff2*adfac
                ad_dmUde(i,j  ,k2)=ad_dmUde(i,j  ,k2)-cff3*adfac
                ad_dmUde(i,j+1,k1)=ad_dmUde(i,j+1,k1)-cff4*adfac
                ad_UFse(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
!>   &                  tl_dZde_p(i,j+1,k1)
!>
                ad_dZde_p(i,j+1,k1)=ad_dZde_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j+1,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
!>   &                  tl_dZde_p(i,j  ,k2)
!>
                ad_dZde_p(i,j  ,k2)=ad_dZde_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j  ,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
!>   &                  tl_dZde_p(i,j+1,k2)
!>
                ad_dZde_p(i,j+1,k2)=ad_dZde_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
!>   &                  tl_dZde_p(i,j  ,k1)
!>
                ad_dZde_p(i,j  ,k1)=ad_dZde_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+     &
!>   &                           cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+     &
!>   &                           cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+     &
!>   &                           cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+             &
     &                   cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+             &
     &                   cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+             &
     &                   cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))*            &
     &                  ad_UFsx(i,j,k2)
#endif
!>              tl_UFsx(i,j,k2)=fac1*                                   &
!>   &                          (tl_cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+  &
!>   &                           tl_cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+  &
!>   &                           tl_cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+  &
!>   &                           tl_cff4*(cff4*dnUdz-dnUdx(i  ,j,k1))+  &
!>   &                           cff1*(tl_cff1*dnUdz+cff1*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i-1,j,k1))+             &
!>   &                           cff2*(tl_cff2*dnUdz+cff2*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i  ,j,k2))+             &
!>   &                           cff3*(tl_cff3*dnUdz+cff3*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i-1,j,k2))+             &
!>   &                           cff4*(tl_cff4*dnUdz+cff4*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i  ,j,k1)))
!>
                cff=2.0_r8*dnUdz
                adfac=fac1*ad_UFsx(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dnUdx(i-1,j,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dnUdx(i  ,j,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dnUdx(i-1,j,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dnUdx(i  ,j,k1))*adfac
                ad_dnUdz=ad_dnUdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dnUdx(i-1,j,k1)=ad_dnUdx(i-1,j,k1)-cff1*adfac
                ad_dnUdx(i  ,j,k2)=ad_dnUdx(i  ,j,k2)-cff2*adfac
                ad_dnUdx(i-1,j,k2)=ad_dnUdx(i-1,j,k2)-cff3*adfac
                ad_dnUdx(i  ,j,k1)=ad_dnUdx(i  ,j,k1)-cff4*adfac
                ad_UFsx(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
!>   &                  tl_dZdx_r(i  ,j,k1)
!>
                ad_dZdx_r(i  ,j,k1)=ad_dZdx_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
!>   &                  tl_dZdx_r(i-1,j,k2)
!>
                ad_dZdx_r(i-1,j,k2)=ad_dZdx_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
!>   &                  tl_dZdx_r(i  ,j,k2)
!>
                ad_dZdx_r(i  ,j,k2)=ad_dZdx_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
!>   &                  tl_dZdx_r(i-1,j,k1)
!>
                ad_dZdx_r(i-1,j,k1)=ad_dZdx_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
!>              tl_dmVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
!>   &                                tl_dVdz(i  ,j+1,k2)+              &
!>   &                                tl_dVdz(i-1,j  ,k2)+              &
!>   &                                tl_dVdz(i  ,j  ,k2))
!>
                adfac=cff*0.25_r8*ad_dmVdz
                ad_dVdz(i-1,j  ,k2)=ad_dVdz(i-1,j  ,k2)+adfac
                ad_dVdz(i  ,j  ,k2)=ad_dVdz(i  ,j  ,k2)+adfac
                ad_dVdz(i-1,j+1,k2)=ad_dVdz(i-1,j+1,k2)+adfac
                ad_dVdz(i  ,j+1,k2)=ad_dVdz(i  ,j+1,k2)+adfac
                ad_dmVdz=0.0_r8
!>              tl_dmUdz=cff*tl_dUdz(i,j,k2)
!>
                ad_dUdz(i,j,k2)=ad_dUdz(i,j,k2)+cff*ad_dmUdz
                ad_dmUdz=0.0_r8
!
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
!>              tl_dnVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
!>   &                                tl_dVdz(i  ,j+1,k2)+              &
!>   &                                tl_dVdz(i-1,j  ,k2)+              &
!>   &                                tl_dVdz(i  ,j  ,k2))
!>
                adfac=cff*0.25_r8*ad_dnVdz
                ad_dVdz(i-1,j  ,k2)=ad_dVdz(i-1,j  ,k2)+adfac
                ad_dVdz(i  ,j  ,k2)=ad_dVdz(i  ,j  ,k2)+adfac
                ad_dVdz(i-1,j+1,k2)=ad_dVdz(i-1,j+1,k2)+adfac
                ad_dVdz(i  ,j+1,k2)=ad_dVdz(i  ,j+1,k2)+adfac
                ad_dnVdz=0.0_r8
!>              tl_dnUdz=cff*tl_dUdz(i,j,k2)
!>
                ad_dUdz(i,j,k2)=ad_dUdz(i,j,k2)+cff*ad_dnUdz
                ad_dnUdz=0.0_r8
#ifdef VISC_3DCOEF
!>              tl_fac2=tl_cff*om_u(i,j)
!>              tl_fac1=tl_cff*on_u(i,j)
!>
                ad_cff=ad_cff+                                          &
     &                 on_u(i,j)*ad_fac1+om_u(i,j)*ad_fac2
                ad_fac1=0.0_r8
                ad_fac2=0.0_r8
# ifdef UV_U3ADV_SPLIT
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_Uvis3d_r(i-1,j,k  )+tl_Uvis3d_r(i,j,k  )+    &
!>   &                  tl_Uvis3d_r(i-1,j,k+1)+tl_Uvis3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_Uvis3d_r(i-1,j,k  )=ad_Uvis3d_r(i-1,j,k  )+adfac
                ad_Uvis3d_r(i  ,j,k  )=ad_Uvis3d_r(i  ,j,k  )+adfac
                ad_Uvis3d_r(i-1,j,k+1)=ad_Uvis3d_r(i-1,j,k+1)+adfac
                ad_Uvis3d_r(i  ,j,k+1)=ad_Uvis3d_r(i  ,j,k+1)+adfac
                ad_cff=0.0_r8
# else
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_visc3d_r(i-1,j,k  )+tl_visc3d_r(i,j,k  )+    &
!>   &                  tl_visc3d_r(i-1,j,k+1)+tl_visc3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_visc3d_r(i-1,j,k  )=ad_visc3d_r(i-1,j,k  )+adfac
                ad_visc3d_r(i  ,j,k  )=ad_visc3d_r(i  ,j,k  )+adfac
                ad_visc3d_r(i-1,j,k+1)=ad_visc3d_r(i-1,j,k+1)+adfac
                ad_visc3d_r(i  ,j,k+1)=ad_visc3d_r(i  ,j,k+1)+adfac
                ad_cff=0.0_r8
# endif
#endif
              END DO
            END DO
          END IF
!
!  Compute adjoint components of the rotated viscous flux (m5/s2) along
!  geopotential surfaces in the XI- and ETA-directions.
!
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
#ifdef VISC_3DCOEF
              cff=0.25_r8*                                              &
     &            (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                           &
     &             Hz(i-1,j-1,k)+Hz(i,j-1,k))*                          &
     &             (on_p(i,j)*(dnVdx(i,j,k1)-                           &
     &                         0.5_r8*pn_p*                             &
     &                         (cff1*(dVdz(i-1,j,k1)+                   &
     &                                dVdz(i  ,j,k2))+                  &
     &                          cff2*(dVdz(i-1,j,k2)+                   &
     &                                dVdz(i  ,j,k1))))+                &
     &              om_p(i,j)*(dmUde(i,j,k1)-                           &
     &                         0.5_r8*pm_p*                             &
     &                         (cff3*(dUdz(i,j-1,k1)+                   &
     &                                dUdz(i,j  ,k2))+                  &
     &                          cff4*(dUdz(i,j-1,k2)+                   &
     &                                dUdz(i,j  ,k1)))))
# ifdef MASKING
              cff=cff*pmask(i,j)
# endif
# ifdef UV_U3ADV_SPLIT
              Uvis_p=0.25_r8*                                           &
     &               (Uvis3d_r(i-1,j-1,k)+Uvis3d_r(i-1,j,k)+            &
     &                Uvis3d_r(i  ,j-1,k)+Uvis3d_r(i  ,j,k))
              Vvis_p=0.25_r8*                                           &
     &               (Vvis3d_r(i-1,j-1,k)+Vvis3d_r(i-1,j,k)+            &
     &                Vvis3d_r(i  ,j-1,k)+Vvis3d_r(i  ,j,k))
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*                          &
!>   &                    (tl_Vvis_p*cff+Vvis_p*tl_cff)
!>
              adfac=on_p(i,j)*on_p(i,j)*ad_VFx(i,j)
              ad_cff=ad_cff+Vvis_p*adfac
              ad_Vvis_p=ad_Vvis_p+cff*adfac
              ad_VFx(i,j)=0.0_r8
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*                          &
!>   &                    (tl_Uvis_p*cff+Uvis_p*tl_cff)
!>
              adfac=om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
              ad_cff=ad_cff+Uvis_p*adfac
              ad_Uvis_p=ad_Uvis_p+cff*adfac
              ad_UFe(i,j)=0.0_r8
!>            tl_Vvis_p=0.25_r8*                                        &
!>   &                  (tl_Vvis3d_r(i-1,j-1,k)+tl_Vvis3d_r(i-1,j,k)+   &
!>   &                   tl_Vvis3d_r(i  ,j-1,k)+tl_Vvis3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_Vvis_p
              ad_Vvis3d_r(i-1,j-1,k)=ad_Vvis3d_r(i-1,j-1,k)+adfac
              ad_Vvis3d_r(i-1,j  ,k)=ad_Vvis3d_r(i-1,j  ,k)+adfac
              ad_Vvis3d_r(i  ,j-1,k)=ad_Vvis3d_r(i  ,j-1,k)+adfac
              ad_Vvis3d_r(i  ,j  ,k)=ad_Vvis3d_r(i  ,j  ,k)+adfac
              ad_Vvis_p=0.0_r8
!>            tl_Uvis_p=0.25_r8*                                        &
!>   &                  (tl_Uvis3d_r(i-1,j-1,k)+tl_Uvis3d_r(i-1,j,k)+   &
!>   &                   tl_Uvis3d_r(i  ,j-1,k)+tl_Uvis3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_Uvis_p
              ad_Uvis3d_r(i-1,j-1,k)=ad_Uvis3d_r(i-1,j-1,k)+adfac
              ad_Uvis3d_r(i-1,j  ,k)=ad_Uvis3d_r(i-1,j  ,k)+adfac
              ad_Uvis3d_r(i  ,j-1,k)=ad_Uvis3d_r(i  ,j-1,k)+adfac
              ad_Uvis3d_r(i  ,j  ,k)=ad_Uvis3d_r(i  ,j  ,k)+adfac
              ad_Uvis_p=0.0_r8
# else
              visc_p=0.25_r8*                                           &
     &               (visc3d_r(i-1,j-1,k)+visc3d_r(i-1,j,k)+            &
     &                visc3d_r(i  ,j-1,k)+visc3d_r(i  ,j,k))
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*                          &
!>   &                    (tl_visc_p*cff+visc_p*tl_cff)
!>
              adfac=on_p(i,j)*on_p(i,j)*ad_VFx(i,j)
              ad_cff=ad_cff+visc_p*adfac
              ad_visc_p=ad_visc_p+cff*adfac
              ad_VFx(i,j)=0.0_r8
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*                          &
!>   &                    (tl_visc_p*cff+visc_p*tl_cff)
!>
              adfac=om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
              ad_cff=ad_cff+visc_p*adfac
              ad_visc_p=ad_visc_p+cff*adfac
              ad_UFe(i,j)=0.0_r8
!>            tl_visc_p=0.25_r8*                                        &
!>   &                  (tl_visc3d_r(i-1,j-1,k)+tl_visc3d_r(i-1,j,k)+   &
!>   &                   tl_visc3d_r(i  ,j-1,k)+tl_visc3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_visc_p
              ad_visc3d_r(i-1,j-1,k)=ad_visc3d_r(i-1,j-1,k)+adfac
              ad_visc3d_r(i-1,j  ,k)=ad_visc3d_r(i-1,j  ,k)+adfac
              ad_visc3d_r(i  ,j-1,k)=ad_visc3d_r(i  ,j-1,k)+adfac
              ad_visc3d_r(i  ,j  ,k)=ad_visc3d_r(i  ,j  ,k)+adfac
              ad_visc_p=0.0_r8
# endif
#else
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*visc4_p(i,j)*tl_cff
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*visc4_p(i,j)*tl_cff
!>
              ad_cff=ad_cff+                                            &
     &               on_p(i,j)*on_p(i,j)*visc4_p(i,j)*ad_VFx(i,j)+      &
     &               om_p(i,j)*om_p(i,j)*visc4_p(i,j)*ad_UFe(i,j)
              ad_VFx(i,j)=0.0_r8
              ad_UFe(i,j)=0.0_r8
#endif
#ifdef MASKING
!>            tl_cff=tl_cff*pmask(i,j)
!>
              ad_cff=ad_cff*pmask(i,j)
#endif
!>            tl_cff=0.25_r8*                                           &
!>   &               ((tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+                 &
!>   &                 tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*                &
!>   &                 (on_p(i,j)*(dnVdx(i,j,k1)-                       &
!>   &                             0.5_r8*pn_p*                         &
!>   &                             (cff1*(dVdz(i-1,j,k1)+               &
!>   &                                    dVdz(i  ,j,k2))+              &
!>   &                              cff2*(dVdz(i-1,j,k2)+               &
!>   &                                    dVdz(i  ,j,k1))))+            &
!>   &                  om_p(i,j)*(dmUde(i,j,k1)-                       &
!>   &                             0.5_r8*pm_p*                         &
!>   &                             (cff3*(dUdz(i,j-1,k1)+               &
!>   &                                    dUdz(i,j  ,k2))+              &
!>   &                              cff4*(dUdz(i,j-1,k2)+               &
!>   &                                    dUdz(i,j  ,k1)))))+           &
!>   &                (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                       &
!>   &                 Hz(i-1,j-1,k)+Hz(i,j-1,k))*                      &
!>   &                 (on_p(i,j)*(tl_dnVdx(i,j,k1)-                    &
!>   &                             0.5_r8*pn_p*                         &
!>   &                             (tl_cff1*(dVdz(i-1,j,k1)+            &
!>   &                                       dVdz(i  ,j,k2))+           &
!>   &                              cff1*(tl_dVdz(i-1,j,k1)+            &
!>   &                                    tl_dVdz(i  ,j,k2))+           &
!>   &                              tl_cff2*(dVdz(i-1,j,k2)+            &
!>   &                                       dVdz(i  ,j,k1))+           &
!>   &                              cff2*(tl_dVdz(i-1,j,k2)+            &
!>   &                                    tl_dVdz(i  ,j,k1))))+         &
!>   &                  om_p(i,j)*(tl_dmUde(i,j,k1)-                    &
!>   &                             0.5_r8*pm_p*                         &
!>   &                             (tl_cff3*(dUdz(i,j-1,k1)+            &
!>   &                                       dUdz(i,j  ,k2))+           &
!>   &                              cff3*(tl_dUdz(i,j-1,k1)+            &
!>   &                                    tl_dUdz(i,j  ,k2))+           &
!>   &                              tl_cff4*(dUdz(i,j-1,k2)+            &
!>   &                                       dUdz(i,j  ,k1))+           &
!>   &                              cff4*(tl_dUdz(i,j-1,k2)+            &
!>   &                                    tl_dUdz(i,j  ,k1))))))
!>
              adfac=0.25_r8*ad_cff
              adfac1=adfac*(on_p(i,j)*(dnVdx(i,j,k1)-                   &
     &                                 0.5_r8*pn_p*                     &
     &                                 (cff1*(dVdz(i-1,j,k1)+           &
     &                                        dVdz(i  ,j,k2))+          &
     &                                  cff2*(dVdz(i-1,j,k2)+           &
     &                                        dVdz(i  ,j,k1))))+        &
     &                      om_p(i,j)*(dmUde(i,j,k1)-                   &
     &                                 0.5_r8*pm_p*                     &
     &                                 (cff3*(dUdz(i,j-1,k1)+           &
     &                                        dUdz(i,j  ,k2))+          &
     &                                  cff4*(dUdz(i,j-1,k2)+           &
     &                                        dUdz(i,j  ,k1)))))
              adfac2=adfac*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+                  &
     &                      Hz(i-1,j-1,k)+Hz(i,j-1,k))
              adfac3=adfac2*on_p(i,j)
              adfac4=adfac3*0.5_r8*pn_p
              adfac5=adfac2*om_p(i,j)
              adfac6=adfac5*0.5_r8*pm_p
              ad_Hz(i-1,j-1,k)=ad_Hz(i-1,j-1,k)+adfac1
              ad_Hz(i  ,j-1,k)=ad_Hz(i  ,j-1,k)+adfac1
              ad_Hz(i-1,j  ,k)=ad_Hz(i-1,j  ,k)+adfac1
              ad_Hz(i  ,j  ,k)=ad_Hz(i  ,j  ,k)+adfac1
              ad_dnVdx(i,j,k1)=ad_dnVdx(i,j,k1)+adfac3
              ad_cff1=ad_cff1-                                          &
     &                (dVdz(i-1,j,k1)+dVdz(i  ,j,k2))*adfac4
              ad_cff2=ad_cff2-                                          &
     &                (dVdz(i-1,j,k2)+dVdz(i  ,j,k1))*adfac4
              ad_dVdz(i-1,j,k1)=ad_dVdz(i-1,j,k1)-cff1*adfac4
              ad_dVdz(i-1,j,k2)=ad_dVdz(i-1,j,k2)-cff2*adfac4
              ad_dVdz(i  ,j,k1)=ad_dVdz(i  ,j,k1)-cff2*adfac4
              ad_dVdz(i  ,j,k2)=ad_dVdz(i  ,j,k2)-cff1*adfac4
              ad_dmUde(i,j,k1)=ad_dmUde(i,j,k1)+adfac5
              ad_cff3=ad_cff3-                                          &
     &                (dUdz(i,j-1,k1)+dUdz(i,j  ,k2))*adfac6
              ad_cff4=ad_cff4-                                          &
     &                (dUdz(i,j-1,k2)+dUdz(i,j  ,k1))*adfac6
              ad_dUdz(i,j-1,k1)=ad_dUdz(i,j-1,k1)-cff3*adfac6
              ad_dUdz(i,j-1,k2)=ad_dUdz(i,j-1,k2)-cff4*adfac6
              ad_dUdz(i,j  ,k1)=ad_dUdz(i,j  ,k1)-cff4*adfac6
              ad_dUdz(i,j  ,k2)=ad_dUdz(i,j  ,k2)-cff3*adfac6
              ad_cff=0.0_r8
!>            tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j,k1)))*            &
!>   &                tl_dZde_p(i,j,k1)
!>            tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j,k1)))*            &
!>   &                tl_dZde_p(i,j,k1)
!>
              ad_dZde_p(i,j,k1)=ad_dZde_p(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZde_p(i,j,k1)))*         &
     &                          ad_cff4+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZde_p(i,j,k1)))*         &
                                ad_cff3
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j,k1)))*            &
!>   &                tl_dZdx_p(i,j,k1)
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j,k1)))*            &
!>   &                tl_dZdx_p(i,j,k1)
!>
              ad_dZdx_p(i,j,k1)=ad_dZdx_p(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZdx_p(i,j,k1)))*         &
     &                          ad_cff2+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZdx_p(i,j,k1)))*         &
                                ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO

          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
              cff3=MIN(dZde_r(i,j,k1),0.0_r8)
              cff4=MAX(dZde_r(i,j,k1),0.0_r8)
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
# ifdef MASKING
              cff=cff*rmask(i,j)
# endif
# ifdef UV_U3ADV_SPLIT
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*                          &
!>   &                    (tl_Vvis3d_r(i,j,k)*cff+                      &
!>   &                     Vvis3d_r(i,j,k)*tl_cff)
!>
              adfac=om_r(i,j)*om_r(i,j)*ad_VFe(i,j)
              ad_cff=ad_cff+Vvis3d_r(i,j,k)*adfac
              ad_Vvis3d_r(i,j,k)=ad_Vvis3d_r(i,j,k)+cff*adfac
              ad_VFe(i,j)=0.0_r8
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*                          &
!>   &                    (tl_Uvis3d_r(i,j,k)*cff+                      &
!>   &                     Uvis3d_r(i,j,k)*tl_cff)
!>
              adfac=on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
              ad_cff=ad_cff+Uvis3d_r(i,j,k)*adfac
              ad_Uvis3d_r(i,j,k)=ad_Uvis3d_r(i,j,k)+cff*adfac
              ad_UFx(i,j)=0.0_r8
# else
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*                          &
!>   &                    (tl_visc3d_r(i,j,k)*cff+                      &
!>   &                     visc3d_r(i,j,k)*tl_cff)
!>
              adfac=om_r(i,j)*om_r(i,j)*ad_VFe(i,j)
              ad_cff=ad_cff+visc3d_r(i,j,k)*adfac
              ad_visc3d_r(i,j,k)=ad_visc3d_r(i,j,k)+cff*adfac
              ad_VFe(i,j)=0.0_r8
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*                          &
!>   &                    (tl_visc3d_r(i,j,k)*cff+                      &
!>   &                     visc3d_r(i,j,k)*tl_cff)
!>
              adfac=on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
              ad_cff=ad_cff+visc3d_r(i,j,k)*adfac
              ad_visc3d_r(i,j,k)=ad_visc3d_r(i,j,k)+cff*adfac
              ad_UFx(i,j)=0.0_r8
# endif
#else
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*visc4_r(i,j)*tl_cff
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*visc4_r(i,j)*tl_cff
!>
              ad_cff=ad_cff+                                            &
     &               om_r(i,j)*om_r(i,j)*visc4_r(i,j)*ad_VFe(i,j)+      &
     &               on_r(i,j)*on_r(i,j)*visc4_r(i,j)*ad_UFx(i,j)
              ad_VFe(i,j)=0.0_r8
              ad_UFx(i,j)=0.0_r8
#endif
#ifdef MASKING
!>            tl_cff=tl_cff*rmask(i,j)
!>
              ad_cff=ad_cff*rmask(i,j)
#endif
!>            tl_cff=tl_Hz(i,j,k)*                                      &
!>   &               (on_r(i,j)*(dnUdx(i,j,k1)-                         &
!>   &                           0.5_r8*pn(i,j)*                        &
!>   &                           (cff1*(dUdz(i  ,j,k1)+                 &
!>   &                                  dUdz(i+1,j,k2))+                &
!>   &                            cff2*(dUdz(i  ,j,k2)+                 &
!>   &                                  dUdz(i+1,j,k1))))-              &
!>   &                om_r(i,j)*(dmVde(i,j,k1)-                         &
!>   &                           0.5_r8*pm(i,j)*                        &
!>   &                           (cff3*(dVdz(i,j  ,k1)+                 &
!>   &                                  dVdz(i,j+1,k2))+                &
!>   &                            cff4*(dVdz(i,j  ,k2)+                 &
!>   &                                  dVdz(i,j+1,k1)))))+             &
!>   &               Hz(i,j,k)*                                         &
!>   &               (on_r(i,j)*(tl_dnUdx(i,j,k1)-                      &
!>   &                           0.5_r8*pn(i,j)*                        &
!>   &                           (tl_cff1*(dUdz(i  ,j,k1)+              &
!>   &                                     dUdz(i+1,j,k2))+             &
!>   &                            cff1*(tl_dUdz(i  ,j,k1)+              &
!>   &                                  tl_dUdz(i+1,j,k2))+             &
!>   &                            tl_cff2*(dUdz(i  ,j,k2)+              &
!>   &                                     dUdz(i+1,j,k1))+             &
!>   &                            cff2*(tl_dUdz(i  ,j,k2)+              &
!>   &                                  tl_dUdz(i+1,j,k1))))-           &
!>   &                om_r(i,j)*(tl_dmVde(i,j,k1)-                      &
!>   &                           0.5_r8*pm(i,j)*                        &
!>   &                           (tl_cff3*(dVdz(i,j  ,k1)+              &
!>   &                                     dVdz(i,j+1,k2))+             &
!>   &                            cff3*(tl_dVdz(i,j  ,k1)+              &
!>   &                                  tl_dVdz(i,j+1,k2))+             &
!>   &                            tl_cff4*(dVdz(i,j  ,k2)+              &
!>   &                                     dVdz(i,j+1,k1))+             &
!>   &                            cff4*(tl_dVdz(i,j  ,k2)+              &
!>   &                                  tl_dVdz(i,j+1,k1)))))
!>
              adfac1=Hz(i,j,k)*ad_cff
              adfac2=adfac1*on_r(i,j)
              adfac3=adfac2*0.5_r8*pn(i,j)
              adfac4=adfac1*om_r(i,j)
              adfac5=adfac4*0.5_r8*pm(i,j)
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                &
     &                     (on_r(i,j)*(dnUdx(i,j,k1)-                   &
     &                                 0.5_r8*pn(i,j)*                  &
     &                                 (cff1*(dUdz(i  ,j,k1)+           &
     &                                        dUdz(i+1,j,k2))+          &
     &                                  cff2*(dUdz(i  ,j,k2)+           &
     &                                        dUdz(i+1,j,k1))))-        &
     &                      om_r(i,j)*(dmVde(i,j,k1)-                   &
     &                                 0.5_r8*pm(i,j)*                  &
     &                                 (cff3*(dVdz(i,j  ,k1)+           &
     &                                        dVdz(i,j+1,k2))+          &
     &                                  cff4*(dVdz(i,j  ,k2)+           &
     &                                        dVdz(i,j+1,k1)))))*       &
     &                     adfac
              ad_dnUdx(i,j,k1)=ad_dnUdx(i,j,k1)+adfac2
              ad_cff1=ad_cff1-                                          &
     &                (dUdz(i  ,j,k1)+dUdz(i+1,j,k2))*adfac3
              ad_cff2=ad_cff2-                                          &
                      (dUdz(i  ,j,k2)+dUdz(i+1,j,k1))*adfac3
              ad_dUdz(i  ,j,k1)=ad_dUdz(i  ,j,k1)-cff1*adfac3
              ad_dUdz(i  ,j,k2)=ad_dUdz(i  ,j,k2)-cff2*adfac3
              ad_dUdz(i+1,j,k1)=ad_dUdz(i+1,j,k1)-cff2*adfac3
              ad_dUdz(i+1,j,k2)=ad_dUdz(i+1,j,k2)-cff1*adfac3
              ad_dmVde(i,j,k1)=ad_dmVde(i,j,k1)-adfac4
              ad_cff3=ad_cff3+                                          &
     &                (dVdz(i,j  ,k1)+dVdz(i,j+1,k2))*adfac5
              ad_cff4=ad_cff4+                                          &
     &                (dVdz(i,j  ,k2)+dVdz(i,j+1,k1))*adfac5
              ad_dVdz(i,j  ,k1)=ad_dVdz(i,j  ,k1)+cff3*adfac5
              ad_dVdz(i,j  ,k2)=ad_dVdz(i,j  ,k2)+cff4*adfac5
              ad_dVdz(i,j+1,k1)=ad_dVdz(i,j+1,k1)+cff4*adfac5
              ad_dVdz(i,j+1,k2)=ad_dVdz(i,j+1,k2)+cff3*adfac5
              ad_cff=0.0_r8
!>            tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j,k1)))*            &
!>   &                tl_dZde_r(i,j,k1)
!>            tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j,k1)))*            &
!>   &                tl_dZde_r(i,j,k1)
!>
              ad_dZde_r(i,j,k1)=ad_dZde_r(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZde_r(i,j,k1)))*         &
     &                          ad_cff4+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZde_r(i,j,k1)))*         &
     &                          ad_cff3
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j,k1)))*            &
!>   &                tl_dZdx_r(i,j,k1)
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j,k1)))*            &
!>   &                tl_dZdx_r(i,j,k1)
!>
              ad_dZdx_r(i,j,k1)=ad_dZdx_r(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZdx_r(i,j,k1)))*         &
     &                          ad_cff2+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZdx_r(i,j,k1)))*         &
     &                          ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO
        END IF
!
!  Compute momentum horizontal (m^-1 s^-3/2) and vertical (s^-3/2)
!  adjoint gradients.
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            tl_VFse(i,j,k2)=0.0_r8
!>
              ad_VFse(i,j,k2)=0.0_r8
!>            tl_VFsx(i,j,k2)=0.0_r8
!>
              ad_VFsx(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            tl_UFse(i,j,k2)=0.0_r8
!>
              ad_UFse(i,j,k2)=0.0_r8
!>            tl_UFsx(i,j,k2)=0.0_r8
!>
              ad_UFsx(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
!>            tl_dVdz(i,j,k2)=0.0_r8
!>
              ad_dVdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
!>            tl_dUdz(i,j,k2)=0.0_r8
!>
              ad_dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-                       &
     &                            z_r(i,j-1,k  )+                       &
     &                            z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j  ,k  )))
!>            tl_dVdz(i,j,k2)=tl_cff*(LapV(i,j,k+1)-                    &
!>   &                                LapV(i,j,k  ))+                   &
!>   &                        cff*(tl_LapV(i,j,k+1)-                    &
!>   &                             tl_LapV(i,j,k  ))
!>
              adfac=cff*ad_dVdz(i,j,k2)
              ad_LapV(i,j,k  )=ad_LapV(i,j,k  )-adfac
              ad_LapV(i,j,k+1)=ad_LapV(i,j,k+1)+adfac
              ad_cff=ad_cff+(LapV(i,j,k+1)-                             &
     &                       LapV(i,j,k  ))*ad_dVdz(i,j,k2)
              ad_dVdz(i,j,k2)=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(tl_z_r(i,j-1,k+1)-               &
!>   &                                 tl_z_r(i,j-1,k  )+               &
!>   &                                 tl_z_r(i,j  ,k+1)-               &
!>   &                                 tl_z_r(i,j  ,k  )))
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-adfac
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO

          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-                       &
     &                            z_r(i-1,j,k  )+                       &
     &                            z_r(i  ,j,k+1)-                       &
     &                            z_r(i  ,j,k  )))
!>            tl_dUdz(i,j,k2)=tl_cff*(LapU(i,j,k+1)-                    &
!>   &                                LapU(i,j,k  ))+                   &
!>   &                        cff*(tl_LapU(i,j,k+1)-                    &
!>   &                             tl_LapU(i,j,k  ))
!>
              adfac=cff*ad_dUdz(i,j,k2)
              ad_LapU(i,j,k  )=ad_LapU(i,j,k  )-adfac
              ad_LapU(i,j,k+1)=ad_LapU(i,j,k+1)+adfac
              ad_cff=ad_cff+(LapU(i,j,k+1)-                             &
     &                       LapU(i,j,k  ))*ad_dUdz(i,j,k2)
              ad_dUdz(i,j,k2)=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*((tl_z_r(i-1,j,k+1)-              &
!>   &                                  tl_z_r(i-1,j,k  )+              &
!>   &                                  tl_z_r(i  ,j,k+1)-              &
!>   &                                  tl_z_r(i  ,j,k  )))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac
              ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-adfac
              ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+adfac
              ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO
        END IF
        IF (k.lt.N(ng)) THEN
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pn(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
!>            tl_dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*              &
!>   &                              tl_LapV(i,j+1,k+1)-                 &
!>   &                              (pm(i,j-1)+pm(i,j  ))*              &
!>   &                              tl_LapV(i,j  ,k+1))
!>
              adfac=cff*ad_dmVde(i,j,k2)
              ad_LapV(i,j  ,k+1)=ad_LapV(i,j  ,k+1)-                    &
     &                           (pm(i,j-1)+pm(i,j  ))*adfac
              ad_LapV(i,j+1,k+1)=ad_LapV(i,j+1,k+1)+                    &
     &                           (pm(i,j  )+pm(i,j+1))*adfac
              ad_dmVde(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstr,Jend+1
            DO i=IstrU-1,Iend+1
              cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                      &
     &                      pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
!>            tl_dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*            &
!>   &                              tl_LapV(i  ,j,k+1)-                 &
!>   &                              (pn(i-1,j-1)+pn(i-1,j))*            &
!>   &                              tl_LapV(i-1,j,k+1))
!>
              adfac=cff*ad_dnVdx(i,j,k2)
              ad_LapV(i-1,j,k+1)=ad_LapV(i-1,j,k+1)-                    &
     &                           (pn(i-1,j-1)+pn(i-1,j))*adfac
              ad_LapV(i  ,j,k+1)=ad_LapV(i  ,j,k+1)+                    &
     &                           (pn(i  ,j-1)+pn(i  ,j))*adfac
              ad_dnVdx(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                      &
     &                      pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
!>            tl_dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*            &
!>   &                              tl_LapU(i,j  ,k+1)-                 &
!>   &                              (pm(i-1,j-1)+pm(i,j-1))*            &
!>   &                              tl_LapU(i,j-1,k+1))
!>
              adfac=cff*ad_dmUde(i,j,k2)
              ad_LapU(i,j-1,k+1)=ad_LapU(i,j-1,k+1)-                    &
     &                           (pm(i-1,j-1)+pm(i,j-1))*adfac
              ad_LapU(i,j  ,k+1)=ad_LapU(i,j  ,k+1)+                    &
     &                           (pm(i-1,j  )+pm(i,j  ))*adfac
              ad_dmUde(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pm(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
!>            tl_dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*              &
!>   &                              tl_LapU(i+1,j,k+1)-                 &
!>   &                              (pn(i-1,j)+pn(i  ,j))*              &
!>   &                              tl_LapU(i  ,j,k+1))
!>
              adfac=cff*ad_dnUdx(i,j,k2)
              ad_LapU(i  ,j,k+1)=ad_LapU(i  ,j,k+1)-                    &
     &                           (pn(i-1,j)+pn(i  ,j))*adfac
              ad_LapU(i+1,j,k+1)=ad_LapU(i+1,j,k+1)+                    &
     &                           (pn(i  ,j)+pn(i+1,j))*adfac
              ad_dnUdx(i,j,k2)=0.0_r8
            END DO
          END DO
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
!>            tl_dZde_r(i,j,k2)=0.5_r8*(tl_VFe(i,j  )+                  &
!>   &                                  tl_VFe(i,j+1))
!>
              adfac=0.5_r8*ad_dZde_r(i,j,k2)
              ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
              ad_VFe(i,j+1)=ad_VFe(i,j+1)+adfac
              ad_dZde_r(i,j,k2)=0.0_r8
!>            tl_dZdx_r(i,j,k2)=0.5_r8*(tl_UFx(i  ,j)+                  &
!>   &                                  tl_UFx(i+1,j))
!>
              adfac=0.5_r8*ad_dZdx_r(i,j,k2)
              ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
              ad_UFx(i+1,j)=ad_UFx(i+1,j)+adfac
              ad_dZdx_r(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
!>            tl_dZde_p(i,j,k2)=0.5_r8*(tl_VFe(i-1,j)+                  &
!>   &                                  tl_VFe(i  ,j))
!>
              adfac=0.5_r8*ad_dZde_p(i,j,k2)
              ad_VFe(i-1,j)=ad_VFe(i-1,j)+adfac
              ad_VFe(i  ,j)=ad_VFe(i  ,j)+adfac
              ad_dZde_p(i,j,k2)=0.0_r8
!>            tl_dZdx_p(i,j,k2)=0.5_r8*(tl_UFx(i,j-1)+                  &
!>   &                                  tl_UFx(i,j  ))
!>
              adfac=0.5_r8*ad_dZdx_p(i,j,k2)
              ad_UFx(i,j-1)=ad_UFx(i,j-1)+adfac
              ad_UFx(i,j  )=ad_UFx(i,j  )+adfac
              ad_dZdx_p(i,j,k2)=0.0_r8
            END DO
          END DO
!
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
!>            tl_VFe(i,j)=cff*(tl_z_r(i,j  ,k+1)-                       &
!>   &                         tl_z_r(i,j-1,k+1))
!>
              adfac=cff*ad_VFe(i,j)
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_VFe(i,j)=0.0_r8
            END DO
          END DO
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
!>            tl_UFx(i,j)=cff*(tl_z_r(i  ,j,k+1)-                       &
!>   &                         tl_z_r(i-1,j,k+1))
!>
              adfac=cff*ad_UFx(i,j)
              ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-adfac
              ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
              ad_UFx(i,j)=0.0_r8
            END DO
          END DO
        END IF
!
!  Compute new storage recursive indices.
!
        kt=k2
        k2=k1
        k1=kt
      END DO K_LOOP2
!
!  Apply boundary conditions (closed or gradient; except periodic)
!  to the first harmonic operator.
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=1,N(ng)
!>          tl_LapV(Iend+1,Jend+1,k)=0.5_r8*                            &
!>   &                               (tl_LapV(Iend  ,Jend+1,k)+         &
!>   &                                tl_LapV(Iend+1,Jend  ,k))
!>
            adfac=0.5_r8*ad_LapV(Iend+1,Jend+1,k)
            ad_LapV(Iend+1,Jend  ,k)=ad_LapV(Iend+1,Jend  ,k)+adfac
            ad_LapV(Iend  ,Jend+1,k)=ad_LapV(Iend  ,Jend+1,k)+adfac
            ad_LapV(Iend+1,Jend+1,k)=0.0_r8
!>          tl_LapU(Iend+1,Jend+1,k)=0.5_r8*                            &
!>   &                               (tl_LapU(Iend  ,Jend+1,k)+         &
!>   &                                tl_LapU(Iend+1,Jend  ,k))
!>
            adfac=0.5_r8*ad_LapU(Iend+1,Jend+1,k)
            ad_LapU(Iend+1,Jend  ,k)=ad_LapU(Iend+1,Jend  ,k)+adfac
            ad_LapU(Iend  ,Jend+1,k)=ad_LapU(Iend  ,Jend+1,k)+adfac
            ad_LapU(Iend+1,Jend+1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=1,N(ng)
!>          tl_LapV(Istr-1,Jend+1,k)=0.5_r8*                            &
!>   &                               (tl_LapV(Istr  ,Jend+1,k)+         &
!>   &                                tl_LapV(Istr-1,Jend  ,k))
!>
            adfac=0.5_r8*ad_LapV(Istr-1,Jend+1,k)
            ad_LapV(Istr-1,Jend  ,k)=ad_LapV(Istr-1,Jend  ,k)+adfac
            ad_LapV(Istr  ,Jend+1,k)=ad_LapV(Istr  ,Jend+1,k)+adfac
            ad_LapV(Istr-1,Jend+1,k)=0.0_r8
!>          tl_LapU(Istr  ,Jend+1,k)=0.5_r8*                            &
!>   &                               (tl_LapU(Istr+1,Jend+1,k)+         &
!>   &                                tl_LapU(Istr  ,Jend  ,k))
!>
            adfac=0.5_r8*ad_LapU(Istr,Jend+1,k)
            ad_LapU(Istr  ,Jend  ,k)=ad_LapU(Istr  ,Jend  ,k)+adfac
            ad_LapU(Istr+1,Jend+1,k)=ad_LapU(Istr+1,Jend+1,k)+adfac
            ad_LapU(Istr  ,Jend+1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
!>          tl_LapV(Iend+1,Jstr  ,k)=0.5_r8*                            &
!>   &                               (tl_LapV(Iend  ,Jstr  ,k)+         &
!>   &                                tl_LapV(Iend+1,Jstr+1,k))
!>
            adfac=0.5_r8*ad_LapV(Iend+1,Jstr,k)
            ad_LapV(Iend  ,Jstr  ,k)=ad_LapV(Iend  ,Jstr  ,k)+adfac
            ad_LapV(Iend+1,Jstr+1,k)=ad_LapV(Iend+1,Jstr+1,k)+adfac
            ad_LapV(Iend+1,Jstr  ,k)=0.0_r8
!>          tl_LapU(Iend+1,Jstr-1,k)=0.5_r8*                            &
!>   &                               (tl_LapU(Iend  ,Jstr-1,k)+         &
!>   &                                tl_LapU(Iend+1,Jstr  ,k))
!>
            adfac=0.5_r8*ad_LapU(Iend+1,Jstr-1,k)
            ad_LapU(Iend  ,Jstr-1,k)=ad_LapU(Iend  ,Jstr-1,k)+adfac
            ad_LapU(Iend+1,Jstr  ,k)=ad_LapU(Iend+1,Jstr  ,k)+adfac
            ad_LapU(Iend+1,Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
!>          tl_LapV(Istr-1,Jstr  ,k)=0.5_r8*                            &
!>   &                               (tl_LapV(Istr-1,Jstr+1,k)+         &
!>   &                                tl_LapV(Istr  ,Jstr  ,k))
!>
            adfac=0.5_r8*ad_LapV(Istr-1,Jstr,k)
            ad_LapV(Istr  ,Jstr  ,k)=ad_LapV(Istr  ,Jstr  ,k)+adfac
            ad_LapV(Istr-1,Jstr+1,k)=ad_LapV(Istr-1,Jstr+1,k)+adfac
            ad_LapV(Istr-1,Jstr  ,k)=0.0_r8
!>          tl_LapU(Istr  ,Jstr-1,k)=0.5_r8*                            &
!>   &                               (tl_LapU(Istr+1,Jstr-1,k)+         &
!>   &                                tl_LapU(Istr  ,Jstr  ,k))
!>
            adfac=0.5_r8*ad_LapU(Istr,Jstr-1,k)
            ad_LapU(Istr+1,Jstr-1,k)=ad_LapU(Istr+1,Jstr-1,k)+adfac
            ad_LapU(Istr  ,Jstr  ,k)=ad_LapU(Istr  ,Jstr  ,k)+adfac
            ad_LapU(Istr  ,Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (ad_LBC(inorth,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,Jend+1,k)=0.0_r8
!>
                ad_LapV(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,Jend+1,k)=tl_LapV(i,Jend,k)
!>
                ad_LapV(i,Jend,k)=ad_LapV(i,Jend,k)+                    &
     &                            ad_LapV(i,Jend+1,k)
                ad_LapV(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(inorth,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jend+1,k)=gamma2(ng)*tl_LapU(i,Jend,k)
!>
                ad_LapU(i,Jend,k)=ad_LapU(i,Jend,k)+                    &
     &                            gamma2(ng)*ad_LapU(i,Jend+1,k)
                ad_LapU(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jend+1,k)=0.0_r8
!>
                 ad_LapU(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (ad_LBC(isouth,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,JstrV-1,k)=0.0_r8
!>
                ad_LapV(i,JstrV-1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,JstrV-1,k)=tl_LapV(i,JstrV,k)
!>
                ad_LapV(i,JstrV,k)=ad_LapV(i,JstrV,k)+                  &
     &                             ad_LapV(i,JstrV-1,k)
                ad_LapV(i,JstrV-1,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(isouth,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jstr-1,k)=gamma2(ng)*tl_LapU(i,Jstr,k)
!>
                ad_LapU(i,Jstr,k)=ad_LapU(i,Jstr,k)+                    &
     &                            gamma2(ng)*ad_LapU(i,Jstr-1,k)
                ad_LapU(i,Jstr-1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jstr-1,k)=0.0_r8
!>
                ad_LapU(i,Jstr-1,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (ad_LBC(ieast,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Iend+1,j,k)=gamma2(ng)*tl_LapV(Iend,j,k)
!>
                ad_LapV(Iend,j,k)=ad_LapV(Iend,j,k)+                    &
     &                            gamma2(ng)*ad_LapV(Iend+1,j,k)
                ad_LapV(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Iend+1,j,k)=0.0_r8
!>
                ad_LapV(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(ieast,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
!>              tl_LapU(Iend+1,j,k)=0.0_r8
!>
                ad_LapU(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
!>              tl_LapU(Iend+1,j,k)=tl_LapU(Iend,j,k)
!>
                ad_LapU(Iend,j,k)=ad_LapU(Iend,j,k)+                    &
     &                            ad_LapU(Iend+1,j,k)
                ad_LapU(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (ad_LBC(iwest,isVvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Istr-1,j,k)=gamma2(ng)*tl_LapV(Istr,j,k)
!>
                ad_LapV(Istr,j,k)=ad_LapV(Istr,j,k)+                    &
     &                            gamma2(ng)*ad_LapV(Istr-1,j,k)
                ad_LapV(Istr-1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Istr-1,j,k)=0.0_r8
!>
                ad_LapV(Istr-1,j,k)=0.0_r8
              END DO
            END DO
          END IF
          IF (ad_LBC(iwest,isUvel,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
!>              tl_LapU(IstrU-1,j,k)=0.0_r8
!>
                ad_LapU(IstrU-1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
!>              tl_LapU(IstrU-1,j,k)=tl_LapU(IstrU,j,k)
!>
                ad_LapU(IstrU,j,k)=ad_LapU(IstrU,j,k)+                  &
     &                             ad_LapU(IstrU-1,j,k)
                ad_LapU(IstrU-1,j,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Compute first adjoint harmonic operator.
!-----------------------------------------------------------------------
!
!  Compute adjoint of starting recursive indices k1 and k2.
!
      k1=2
      k2=1
      DO k=0,N(ng)
!!
!!  Note: The following code is equivalent to
!!
!!      kt=k1
!!      k1=k2
!!      k2=kt
!!
!!  We use the adjoint of above code.
!!
        k1=k2
        k2=3-k1
      END DO
!
!  Compute required BASIC STATE fields. Need to look forward in "kk"
!  index.
!
      K_LOOP3: DO k=N(ng),0,-1
        k2b=1
        DO kk=0,k
          k1b=k2b
          k2b=3-k1b
          IF (kk.lt.N(ng)) THEN
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
            DO j=Jstrm2,Jendp2
              DO i=IstrUm2,Iendp2
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
                UFx(i,j)=cff*(z_r(i  ,j,kk+1)-                          &
     &                        z_r(i-1,j,kk+1))
              END DO
            END DO
            DO j=JstrVm2,Jendp2
              DO i=Istrm2,Iendp2
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
                VFe(i,j)=cff*(z_r(i,j  ,kk+1)-                          &
     &                        z_r(i,j-1,kk+1))
              END DO
            END DO
!
            DO j=Jstrm1,Jendp2
              DO i=Istrm1,Iendp2
                dZdx_p(i,j,k2b)=0.5_r8*(UFx(i,j-1)+                     &
     &                                  UFx(i,j  ))
                dZde_p(i,j,k2b)=0.5_r8*(VFe(i-1,j)+                     &
     &                                  VFe(i  ,j))
              END DO
            END DO
            DO j=JstrVm2,Jendp1
              DO i=IstrUm2,Iendp1
                dZdx_r(i,j,k2b)=0.5_r8*(UFx(i  ,j)+                     &
     &                                  UFx(i+1,j))
                dZde_r(i,j,k2b)=0.5_r8*(VFe(i,j  )+                     &
     &                                  VFe(i,j+1))
              END DO
            END DO
!
!  Compute momentum horizontal (1/m/s) and vertical (1/s) gradients.
!
            DO j=JstrVm2,Jendp1
              DO i=IstrUm2,Iendp1
                cff=0.5_r8*pm(i,j)
#ifdef MASKING
                cff=cff*rmask(i,j)
#endif
                dnUdx(i,j,k2b)=cff*((pn(i  ,j)+pn(i+1,j))*              &
     &                              u(i+1,j,kk+1,nrhs)-                 &
     &                              (pn(i-1,j)+pn(i  ,j))*              &
     &                              u(i  ,j,kk+1,nrhs))
              END DO
            END DO

            DO j=Jstrm1,Jendp2
              DO i=Istrm1,Iendp2
                cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                    &
     &                        pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
                cff=cff*pmask(i,j)
#endif
                dmUde(i,j,k2b)=cff*((pm(i-1,j  )+pm(i,j  ))*            &
     &                              u(i,j  ,kk+1,nrhs)-                 &
     &                              (pm(i-1,j-1)+pm(i,j-1))*            &
     &                              u(i,j-1,kk+1,nrhs))
              END DO
            END DO

            DO j=Jstrm1,Jendp2
              DO i=Istrm1,Iendp2
                cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                    &
     &                        pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
                cff=cff*pmask(i,j)
#endif
                dnVdx(i,j,k2b)=cff*((pn(i  ,j-1)+pn(i  ,j))*            &
     &                              v(i  ,j,kk+1,nrhs)-                 &
     &                              (pn(i-1,j-1)+pn(i-1,j))*            &
     &                              v(i-1,j,kk+1,nrhs))
              END DO
            END DO

            DO j=JstrVm2,Jendp1
              DO i=IstrUm2,Iendp1
                cff=0.5_r8*pn(i,j)
#ifdef MASKING
                cff=cff*rmask(i,j)
#endif
                dmVde(i,j,k2b)=cff*((pm(i,j  )+pm(i,j+1))*              &
     &                              v(i,j+1,kk+1,nrhs)-                 &
     &                              (pm(i,j-1)+pm(i,j  ))*              &
     &                              v(i,j  ,kk+1,nrhs))
              END DO
            END DO
          END IF

          IF ((kk.eq.0).or.(kk.eq.N(ng))) THEN
            DO j=Jstrm2,Jendp2
              DO i=IstrUm2,Iendp2
                dUdz(i,j,k2b)=0.0_r8
              END DO
            END DO
            DO j=JstrVm2,Jendp2
              DO i=Istrm2,Iendp2
                dVdz(i,j,k2b)=0.0_r8
              END DO
            END DO

            DO j=Jstrm1,Jendp1
              DO i=IstrUm1,Iendp1
                UFsx(i,j,k2b)=0.0_r8
                UFse(i,j,k2b)=0.0_r8
              END DO
            END DO
            DO j=JstrVm1,Jendp1
              DO i=Istrm1,Iendp1
                VFsx(i,j,k2b)=0.0_r8
                VFse(i,j,k2b)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstrm2,Jendp2
              DO i=IstrUm2,Iendp2
                cff=1.0_r8/(0.5_r8*(z_r(i-1,j,kk+1)-                    &
     &                              z_r(i-1,j,kk  )+                    &
     &                              z_r(i  ,j,kk+1)-                    &
     &                              z_r(i  ,j,kk  )))
                dUdz(i,j,k2b)=cff*(u(i,j,kk+1,nrhs)-                    &
     &                             u(i,j,kk  ,nrhs))
              END DO
            END DO

            DO j=JstrVm2,Jendp2
              DO i=Istrm2,Iendp2
                cff=1.0_r8/(0.5_r8*(z_r(i,j-1,kk+1)-                    &
     &                              z_r(i,j-1,kk  )+                    &
     &                              z_r(i,j  ,kk+1)-                    &
     &                              z_r(i,j  ,kk  )))
                dVdz(i,j,k2b)=cff*(v(i,j,kk+1,nrhs)-                    &
     &                             v(i,j,kk  ,nrhs))
              END DO
            END DO
          END IF
!
!  Compute BASIC STATE components of the rotated viscous flux
!  (m^4 s-^3/2) along geopotential surfaces in the XI- and
!  ETA-directions.
!
          IF (kk.gt.0) THEN
            DO j=JstrVm2,Jendp1
              DO i=IstrUm2,Iendp1
                cff1=MIN(dZdx_r(i,j,k1b),0.0_r8)
                cff2=MAX(dZdx_r(i,j,k1b),0.0_r8)
                cff3=MIN(dZde_r(i,j,k1b),0.0_r8)
                cff4=MAX(dZde_r(i,j,k1b),0.0_r8)
                cff=on_r(i,j)*(dnUdx(i,j,k1b)-                          &
     &                         0.5_r8*pn(i,j)*                          &
     &                         (cff1*(dUdz(i  ,j,k1b)+                  &
     &                                dUdz(i+1,j,k2b))+                 &
     &                          cff2*(dUdz(i  ,j,k2b)+                  &
     &                                dUdz(i+1,j,k1b))))-               &
     &              om_r(i,j)*(dmVde(i,j,k1b)-                          &
     &                         0.5_r8*pm(i,j)*                          &
     &                         (cff3*(dVdz(i,j  ,k1b)+                  &
     &                                dVdz(i,j+1,k2b))+                 &
     &                          cff4*(dVdz(i,j  ,k2b)+                  &
     &                                dVdz(i,j+1,k1b))))
#ifdef MASKING
                cff=cff*rmask(i,j)
#endif
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                UFx(i,j)=on_r(i,j)*on_r(i,j)*Uvis3d_r(i,j,kk)*cff
                VFe(i,j)=om_r(i,j)*om_r(i,j)*Vvis3d_r(i,j,kk)*cff
# else
                UFx(i,j)=on_r(i,j)*on_r(i,j)*visc3d_r(i,j,kk)*cff
                VFe(i,j)=om_r(i,j)*om_r(i,j)*visc3d_r(i,j,kk)*cff
# endif
#else
                UFx(i,j)=on_r(i,j)*on_r(i,j)*visc4_r(i,j)*cff
                VFe(i,j)=om_r(i,j)*om_r(i,j)*visc4_r(i,j)*cff
#endif
              END DO
            END DO

            DO j=Jstrm1,Jendp2
              DO i=Istrm1,Iendp2
                pm_p=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                    &
     &                        pm(i  ,j-1)+pm(i  ,j))
                pn_p=0.25_r8*(pn(i-1,j-1)+pn(i-1,j)+                    &
     &                        pn(i  ,j-1)+pn(i  ,j))
                cff1=MIN(dZdx_p(i,j,k1b),0.0_r8)
                cff2=MAX(dZdx_p(i,j,k1b),0.0_r8)
                cff3=MIN(dZde_p(i,j,k1b),0.0_r8)
                cff4=MAX(dZde_p(i,j,k1b),0.0_r8)
                cff=on_p(i,j)*(dnVdx(i,j,k1b)-                          &
     &                         0.5_r8*pn_p*                             &
     &                         (cff1*(dVdz(i-1,j,k1b)+                  &
     &                                dVdz(i  ,j,k2b))+                 &
     &                          cff2*(dVdz(i-1,j,k2b)+                  &
     &                                dVdz(i  ,j,k1b))))+               &
     &              om_p(i,j)*(dmUde(i,j,k1b)-                          &
     &                         0.5_r8*pm_p*                             &
     &                         (cff3*(dUdz(i,j-1,k1b)+                  &
     &                                dUdz(i,j  ,k2b))+                 &
     &                          cff4*(dUdz(i,j-1,k2b)+                  &
     &                                dUdz(i,j  ,k1b))))
#ifdef MASKING
                cff=cff*pmask(i,j)
#endif
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                Uvis_p=0.25_r8*                                         &
     &                 (Uvis3d_r(i-1,j-1,kk)+Uvis3d_r(i-1,j,kk)+        &
     &                  Uvis3d_r(i  ,j-1,kk)+Uvis3d_r(i  ,j,kk))
                Vvis_p=0.25_r8*                                         &
     &                 (Vvis3d_r(i-1,j-1,kk)+Vvis3d_r(i-1,j,kk)+        &
     &                  Vvis3d_r(i  ,j-1,kk)+Vvis3d_r(i  ,j,kk))
                UFe(i,j)=om_p(i,j)*om_p(i,j)*Uvis_p*cff
                VFx(i,j)=on_p(i,j)*on_p(i,j)*Vvis_p*cff
# else
                visc_p=0.25_r8*                                         &
     &                 (visc3d_r(i-1,j-1,kk)+visc3d_r(i-1,j,kk)+        &
     &                  visc3d_r(i  ,j-1,kk)+visc3d_r(i  ,j,kk))
                UFe(i,j)=om_p(i,j)*om_p(i,j)*visc_p*cff
                VFx(i,j)=on_p(i,j)*on_p(i,j)*visc_p*cff
# endif
#else
                UFe(i,j)=om_p(i,j)*om_p(i,j)*visc4_p(i,j)*cff
                VFx(i,j)=on_p(i,j)*on_p(i,j)*visc4_p(i,j)*cff
#endif
              END DO
            END DO
!
!  Compute BASIC STATE vertical flux (m^2 s^-3/2) due to sloping
!  terrain-following surfaces.
!
            IF (kk.lt.N(ng)) THEN
              DO j=JstrVm2,Jendp1
                DO i=IstrUm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                  cff=0.125_r8*                                         &
     &                (Uvis3d_r(i-1,j,kk  )+Uvis3d_r(i,j,kk  )+         &
     &                 Uvis3d_r(i-1,j,kk+1)+Uvis3d_r(i,j,kk+1))
# else
                  cff=0.125_r8*                                         &
     &                (visc3d_r(i-1,j,kk  )+visc3d_r(i,j,kk  )+         &
     &                 visc3d_r(i-1,j,kk+1)+visc3d_r(i,j,kk+1))
# endif
                  fac1=cff*on_u(i,j)
                  fac2=cff*om_u(i,j)
#else
                  cff=0.25_r8*(visc4_r(i-1,j)+visc4_r(i,j))
                  fac1=cff*on_u(i,j)
                  fac2=cff*om_u(i,j)
#endif
                  cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                  dnUdz=cff*dUdz(i,j,k2b)
                  dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2b)+                 &
     &                               dVdz(i  ,j+1,k2b)+                 &
     &                               dVdz(i-1,j  ,k2b)+                 &
     &                               dVdz(i  ,j  ,k2b))
                  cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                  dmUdz=cff*dUdz(i,j,k2b)
                  dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2b)+                 &
     &                               dVdz(i  ,j+1,k2b)+                 &
     &                               dVdz(i-1,j  ,k2b)+                 &
     &                               dVdz(i  ,j  ,k2b))

                  cff1=MIN(dZdx_r(i-1,j,k1b),0.0_r8)
                  cff2=MIN(dZdx_r(i  ,j,k2b),0.0_r8)
                  cff3=MAX(dZdx_r(i-1,j,k2b),0.0_r8)
                  cff4=MAX(dZdx_r(i  ,j,k1b),0.0_r8)
                  UFsx(i,j,k2b)=fac1*                                   &
     &                          (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1b))+    &
     &                           cff2*(cff2*dnUdz-dnUdx(i  ,j,k2b))+    &
     &                           cff3*(cff3*dnUdz-dnUdx(i-1,j,k2b))+    &
     &                           cff4*(cff4*dnUdz-dnUdx(i  ,j,k1b)))

                  cff1=MIN(dZde_p(i,j  ,k1b),0.0_r8)
                  cff2=MIN(dZde_p(i,j+1,k2b),0.0_r8)
                  cff3=MAX(dZde_p(i,j  ,k2b),0.0_r8)
                  cff4=MAX(dZde_p(i,j+1,k1b),0.0_r8)
                  UFse(i,j,k2b)=fac2*                                   &
     &                          (cff1*(cff1*dmUdz-dmUde(i,j  ,k1b))+    &
     &                           cff2*(cff2*dmUdz-dmUde(i,j+1,k2b))+    &
     &                           cff3*(cff3*dmUdz-dmUde(i,j  ,k2b))+    &
     &                           cff4*(cff4*dmUdz-dmUde(i,j+1,k1b)))

                  cff1=MIN(dZde_p(i,j  ,k1b),0.0_r8)
                  cff2=MIN(dZde_p(i,j+1,k2b),0.0_r8)
                  cff3=MAX(dZde_p(i,j  ,k2b),0.0_r8)
                  cff4=MAX(dZde_p(i,j+1,k1b),0.0_r8)
                  cff5=MIN(dZdx_p(i,j  ,k1b),0.0_r8)
                  cff6=MIN(dZdx_p(i,j+1,k2b),0.0_r8)
                  cff7=MAX(dZdx_p(i,j  ,k2b),0.0_r8)
                  cff8=MAX(dZdx_p(i,j+1,k1b),0.0_r8)
                  UFsx(i,j,k2b)=UFsx(i,j,k2b)+                          &
     &                          fac1*                                   &
     &                          (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1b))+    &
     &                           cff2*(cff6*dnVdz-dnVdx(i,j+1,k2b))+    &
     &                           cff3*(cff7*dnVdz-dnVdx(i,j  ,k2b))+    &
     &                           cff4*(cff8*dnVdz-dnVdx(i,j+1,k1b)))

                  cff1=MIN(dZdx_r(i-1,j,k1b),0.0_r8)
                  cff2=MIN(dZdx_r(i  ,j,k2b),0.0_r8)
                  cff3=MAX(dZdx_r(i-1,j,k2b),0.0_r8)
                  cff4=MAX(dZdx_r(i  ,j,k1b),0.0_r8)
                  cff5=MIN(dZde_r(i-1,j,k1b),0.0_r8)
                  cff6=MIN(dZde_r(i  ,j,k2b),0.0_r8)
                  cff7=MAX(dZde_r(i-1,j,k2b),0.0_r8)
                  cff8=MAX(dZde_r(i  ,j,k1b),0.0_r8)
                  UFse(i,j,k2b)=UFse(i,j,k2b)-                          &
     &                          fac2*                                   &
     &                          (cff1*(cff5*dmVdz-dmVde(i-1,j,k1b))+    &
     &                           cff2*(cff6*dmVdz-dmVde(i  ,j,k2b))+    &
     &                           cff3*(cff7*dmVdz-dmVde(i-1,j,k2b))+    &
     &                           cff4*(cff8*dmVdz-dmVde(i  ,j,k1b)))
                END DO
              END DO
!
              DO j=JstrVm1,Jendp1
                DO i=Istrm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                  cff=0.125_r8*                                         &
     &                (Vvis3d_r(i,j-1,kk  )+Vvis3d_r(i,j,kk  )+         &
     &                 Vvis3d_r(i,j-1,kk+1)+Vvis3d_r(i,j,kk+1))
# else
                  cff=0.125_r8*                                         &
     &                (visc3d_r(i,j-1,kk  )+visc3d_r(i,j,kk  )+         &
     &                 visc3d_r(i,j-1,kk+1)+visc3d_r(i,j,kk+1))
# endif
                  fac1=cff*on_v(i,j)
                  fac2=cff*om_v(i,j)
#else
                  cff=0.25_r8*(visc4_r(i,j-1)+visc4_r(i,j))
                  fac1=cff*on_v(i,j)
                  fac2=cff*om_v(i,j)
#endif
                  cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                  dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2b)+                 &
     &                               dUdz(i+1,j  ,k2b)+                 &
     &                               dUdz(i  ,j-1,k2b)+                 &
     &                               dUdz(i+1,j-1,k2b))
                  dnVdz=cff*dVdz(i,j,k2b)
                  cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                  dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2b)+                 &
     &                               dUdz(i+1,j  ,k2b)+                 &
     &                               dUdz(i  ,j-1,k2b)+                 &
     &                               dUdz(i+1,j-1,k2b))
                  dmVdz=cff*dVdz(i,j,k2b)

                  cff1=MIN(dZdx_p(i  ,j,k1b),0.0_r8)
                  cff2=MIN(dZdx_p(i+1,j,k2b),0.0_r8)
                  cff3=MAX(dZdx_p(i  ,j,k2b),0.0_r8)
                  cff4=MAX(dZdx_p(i+1,j,k1b),0.0_r8)
                  VFsx(i,j,k2b)=fac1*                                   &
     &                          (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1b))+    &
     &                           cff2*(cff2*dnVdz-dnVdx(i+1,j,k2b))+    &
     &                           cff3*(cff3*dnVdz-dnVdx(i  ,j,k2b))+    &
     &                           cff4*(cff4*dnVdz-dnVdx(i+1,j,k1b)))

                  cff1=MIN(dZde_r(i,j-1,k1b),0.0_r8)
                  cff2=MIN(dZde_r(i,j  ,k2b),0.0_r8)
                  cff3=MAX(dZde_r(i,j-1,k2b),0.0_r8)
                  cff4=MAX(dZde_r(i,j  ,k1b),0.0_r8)
                  VFse(i,j,k2b)=fac2*                                   &
     &                          (cff1*(cff1*dmVdz-dmVde(i,j-1,k1b))+    &
     &                           cff2*(cff2*dmVdz-dmVde(i,j  ,k2b))+    &
     &                           cff3*(cff3*dmVdz-dmVde(i,j-1,k2b))+    &
     &                           cff4*(cff4*dmVdz-dmVde(i,j  ,k1b)))

                  cff1=MIN(dZde_r(i,j-1,k1b),0.0_r8)
                  cff2=MIN(dZde_r(i,j  ,k2b),0.0_r8)
                  cff3=MAX(dZde_r(i,j-1,k2b),0.0_r8)
                  cff4=MAX(dZde_r(i,j  ,k1b),0.0_r8)
                  cff5=MIN(dZdx_r(i,j-1,k1b),0.0_r8)
                  cff6=MIN(dZdx_r(i,j  ,k2b),0.0_r8)
                  cff7=MAX(dZdx_r(i,j-1,k2b),0.0_r8)
                  cff8=MAX(dZdx_r(i,j  ,k1b),0.0_r8)
                  VFsx(i,j,k2b)=VFsx(i,j,k2b)-                          &
     &                          fac1*                                   &
     &                          (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1b))+    &
     &                           cff2*(cff6*dnUdz-dnUdx(i,j  ,k2b))+    &
     &                           cff3*(cff7*dnUdz-dnUdx(i,j-1,k2b))+    &
     &                           cff4*(cff8*dnUdz-dnUdx(i,j  ,k1b)))

                  cff1=MIN(dZdx_p(i  ,j,k1b),0.0_r8)
                  cff2=MIN(dZdx_p(i+1,j,k2b),0.0_r8)
                  cff3=MAX(dZdx_p(i  ,j,k2b),0.0_r8)
                  cff4=MAX(dZdx_p(i+1,j,k1b),0.0_r8)
                  cff5=MIN(dZde_p(i  ,j,k1b),0.0_r8)
                  cff6=MIN(dZde_p(i+1,j,k2b),0.0_r8)
                  cff7=MAX(dZde_p(i  ,j,k2b),0.0_r8)
                  cff8=MAX(dZde_p(i+1,j,k1b),0.0_r8)
                  VFse(i,j,k2b)=VFse(i,j,k2b)+                          &
     &                          fac2*                                   &
     &                          (cff1*(cff5*dmUdz-dmUde(i  ,j,k1b))+    &
     &                           cff2*(cff6*dmUdz-dmUde(i+1,j,k2b))+    &
     &                           cff3*(cff7*dmUdz-dmUde(i  ,j,k2b))+    &
     &                           cff4*(cff8*dmUdz-dmUde(i+1,j,k1b)))
                END DO
              END DO
            END IF
          END IF
        END DO
!
        IF (k.gt.0) THEN
!
!  Compute first adjoint harmonic operator (m s^-3/2).
!
          DO j=JstrVm1,Jendp1
            DO i=Istrm1,Iendp1
              cff=0.125_r8*(pm(i,j)+pm(i,j-1))*                         &
     &                     (pn(i,j)+pn(i,j-1))
              cff1=1.0_r8/(0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k)))
#ifdef MASKING
!>            tl_LapV(i,j,k)=tl_LapV(i,j,k)*vmask(i,j)
!>
              ad_LapV(i,j,k)=ad_LapV(i,j,k)*vmask(i,j)
#endif
!>            tl_LapV(i,j,k)=cff*((pn(i,j-1)+pn(i,j))*                  &
!>   &                            (tl_VFx(i+1,j)-tl_VFx(i,j))-          &
!>   &                            (pm(i,j-1)+pm(i,j))*                  &
!>   &                            (tl_VFe(i,j)-tl_VFe(i,j-1)))+         &
!>   &                       tl_cff1*((VFsx(i,j,k2)+VFse(i,j,k2))-      &
!>   &                                (VFsx(i,j,k1)+VFse(i,j,k1)))+     &
!>   &                       cff1*((tl_VFsx(i,j,k2)+tl_VFse(i,j,k2))-   &
!>   &                             (tl_VFsx(i,j,k1)+tl_VFse(i,j,k1)))
!>
              adfac=cff1*ad_LapV(i,j,k)
              adfac1=cff*ad_LapV(i,j,k)
              adfac2=adfac1*(pm(i,j-1)+pm(i,j))
              adfac3=adfac1*(pn(i,j-1)+pn(i,j))
              ad_VFsx(i,j,k1)=ad_VFsx(i,j,k1)-adfac
              ad_VFse(i,j,k1)=ad_VFse(i,j,k1)-adfac
              ad_VFsx(i,j,k2)=ad_VFsx(i,j,k2)+adfac
              ad_VFse(i,j,k2)=ad_VFse(i,j,k2)+adfac
              ad_cff1=ad_cff1+                                          &
     &                ((VFsx(i,j,k2)+VFse(i,j,k2))-                     &
     &                 (VFsx(i,j,k1)+VFse(i,j,k1)))*ad_LapV(i,j,k)
              ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac2
              ad_VFe(i,j  )=ad_VFe(i,j  )-adfac2
              ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac3
              ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac3
              ad_LapV(i,j,k)=0.0_r8
!>            tl_cff1=-cff1*cff1*(0.5_r8*(tl_Hz(i,j-1,k)+tl_Hz(i,j,k)))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
              ad_cff1=0.0_r8
            END DO
          END DO

          DO j=Jstrm1,Jendp1
            DO i=IstrUm1,Iendp1
              cff=0.125_r8*(pm(i-1,j)+pm(i,j))*                         &
     &                     (pn(i-1,j)+pn(i,j))
              cff1=1.0_r8/(0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k)))
#ifdef MASKING
!>            tl_LapU(i,j,k)=tl_LapU(i,j,k)*umask(i,j)
!>
              ad_LapU(i,j,k)=ad_LapU(i,j,k)*umask(i,j)
#endif
!>            tl_LapU(i,j,k)=cff*((pn(i-1,j)+pn(i,j))*                  &
!>                                (tl_UFx(i,j)-tl_UFx(i-1,j))+          &
!>   &                            (pm(i-1,j)+pm(i,j))*                  &
!>   &                            (tl_UFe(i,j+1)-tl_UFe(i,j)))+         &
!>   &                       tl_cff1*((UFsx(i,j,k2)+UFse(i,j,k2))-      &
!>   &                                (UFsx(i,j,k1)+UFse(i,j,k1)))+     &
!>   &                       cff1*((tl_UFsx(i,j,k2)+tl_UFse(i,j,k2))-   &
!>   &                             (tl_UFsx(i,j,k1)+tl_UFse(i,j,k1)))
!>
              adfac=cff1*ad_LapU(i,j,k)
              adfac1=cff*ad_LapU(i,j,k)
              adfac2=adfac1*(pm(i-1,j)+pm(i,j))
              adfac3=adfac1*(pn(i-1,j)+pn(i,j))
              ad_UFsx(i,j,k1)=ad_UFsx(i,j,k1)-adfac
              ad_UFse(i,j,k1)=ad_UFse(i,j,k1)-adfac
              ad_UFsx(i,j,k2)=ad_UFsx(i,j,k2)+adfac
              ad_UFse(i,j,k2)=ad_UFse(i,j,k2)+adfac
              ad_cff1=ad_cff1+                                          &
     &                ((UFsx(i,j,k2)+UFse(i,j,k2))-                     &
     &                 (UFsx(i,j,k1)+UFse(i,j,k1)))*ad_LapU(i,j,k)
              ad_UFe(i,j  )=ad_UFe(i,j  )-adfac2
              ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac2
              ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac3
              ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac3
              ad_LapU(i,j,k)=0.0_r8
!>            tl_cff1=-cff1*cff1*(0.5_r8*(tl_Hz(i-1,j,k)+tl_Hz(i,j,k)))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
              ad_cff1=0.0_r8
            END DO
          END DO
!
!  Compute vertical flux (m^2 s^-3/2) due to sloping terrain-following
!  surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=JstrVm1,Jendp1
              DO i=Istrm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Vvis3d_r(i,j-1,k  )+Vvis3d_r(i,j,k  )+             &
     &               Vvis3d_r(i,j-1,k+1)+Vvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i,j-1,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i,j-1,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#else
                cff=0.25_r8*(visc4_r(i,j-1)+visc4_r(i,j))
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#endif
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)
!
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+     &
!>   &                           cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+     &
!>   &                           cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+     &
!>   &                           cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+             &
     &                   cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+             &
     &                   cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+             &
     &                   cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))*            &
     &                  ad_VFse(i,j,k2)
#endif
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          fac2*                                   &
!>   &                          (tl_cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+  &
!>   &                           tl_cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+  &
!>   &                           tl_cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+  &
!>   &                           tl_cff4*(cff8*dmUdz-dmUde(i+1,j,k1))+  &
!>   &                           cff1*(tl_cff5*dmUdz+cff5*tl_dmUdz-     &
!>   &                                 tl_dmUde(i  ,j,k1))+             &
!>   &                           cff2*(tl_cff6*dmUdz+cff6*tl_dmUdz-     &
!>   &                                 tl_dmUde(i+1,j,k2))+             &
!>   &                           cff3*(tl_cff7*dmUdz+cff7*tl_dmUdz-     &
!>   &                                 tl_dmUde(i  ,j,k2))+             &
!>   &                           cff4*(tl_cff8*dmUdz+cff8*tl_dmUdz-     &
!>   &                                 tl_dmUde(i+1,j,k1)))
!>
                adfac=fac2*ad_VFse(i,j,k2)
                adfac1=adfac*dmUdz
                ad_cff1=ad_cff1+(cff5*dmUdz-dmUde(i  ,j,k1))*adfac
                ad_cff2=ad_cff2+(cff6*dmUdz-dmUde(i+1,j,k2))*adfac
                ad_cff3=ad_cff3+(cff7*dmUdz-dmUde(i  ,j,k2))*adfac
                ad_cff4=ad_cff4+(cff8*dmUdz-dmUde(i+1,j,k1))*adfac
                ad_cff5=ad_cff5+cff1*adfac1
                ad_cff6=ad_cff6+cff2*adfac1
                ad_cff7=ad_cff7+cff3*adfac1
                ad_cff8=ad_cff8+cff4*adfac1
                ad_dmUdz=ad_dmUdz+                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dmUde(i  ,j,k1)=ad_dmUde(i  ,j,k1)-cff1*adfac
                ad_dmUde(i+1,j,k2)=ad_dmUde(i+1,j,k2)-cff2*adfac
                ad_dmUde(i  ,j,k2)=ad_dmUde(i  ,j,k2)-cff3*adfac
                ad_dmUde(i+1,j,k1)=ad_dmUde(i+1,j,k1)-cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_p(i+1,j,k1)))*        &
!>   &                  tl_dZde_p(i+1,j,k1)
!>
                ad_dZde_p(i+1,j,k1)=ad_dZde_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i+1,j,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_p(i  ,j,k2)))*        &
!>   &                  tl_dZde_p(i  ,j,k2)
!>
                ad_dZde_p(i  ,j,k2)=ad_dZde_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i  ,j,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_p(i+1,j,k2)))*        &
!>   &                  tl_dZde_p(i+1,j,k2)
!>
                ad_dZde_p(i+1,j,k2)=ad_dZde_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i+1,j,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_p(i  ,j,k1)))*        &
!>   &                  tl_dZde_p(i  ,j,k1)
!>
                ad_dZde_p(i  ,j,k1)=ad_dZde_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i  ,j,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
!>   &                  tl_dZdx_p(i+1,j,k1)
!>
                ad_dZdx_p(i+1,j,k1)=ad_dZdx_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
!>   &                  tl_dZdx_p(i  ,j,k2)
!>
                ad_dZdx_p(i  ,j,k2)=ad_dZdx_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
!>   &                  tl_dZdx_p(i+1,j,k2)
!>
                ad_dZdx_p(i+1,j,k2)=ad_dZdx_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
!>   &                  tl_dZdx_p(i  ,j,k1)
!>
                ad_dZdx_p(i  ,j,k1)=ad_dZdx_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+     &
!>   &                           cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+     &
!>   &                           cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+     &
!>   &                           cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
!>
                ad_fac1=ad_fac1-                                        &
     &                  (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+             &
     &                   cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+             &
     &                   cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+             &
     &                   cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))*            &
     &                  ad_VFsx(i,j,k2)
#endif
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
!>   &                          fac1*                                   &
!>   &                          (tl_cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+  &
!>   &                           tl_cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+  &
!>   &                           tl_cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+  &
!>   &                           tl_cff4*(cff8*dnUdz-dnUdx(i,j  ,k1))+  &
!>   &                           cff1*(tl_cff5*dnUdz+cff5*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j-1,k1))+             &
!>   &                           cff2*(tl_cff6*dnUdz+cff6*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j  ,k2))+             &
!>   &                           cff3*(tl_cff7*dnUdz+cff7*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j-1,k2))+             &
!>   &                           cff4*(tl_cff8*dnUdz+cff8*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i,j  ,k1)))
!>
                adfac=fac1*ad_VFsx(i,j,k2)
                adfac1=adfac*dnUdz
                ad_cff1=ad_cff1-(cff5*dnUdz-dnUdx(i,j-1,k1))*adfac
                ad_cff2=ad_cff2-(cff6*dnUdz-dnUdx(i,j  ,k2))*adfac
                ad_cff3=ad_cff3-(cff7*dnUdz-dnUdx(i,j-1,k2))*adfac
                ad_cff4=ad_cff4-(cff8*dnUdz-dnUdx(i,j  ,k1))*adfac
                ad_cff5=ad_cff5-cff1*adfac1
                ad_cff6=ad_cff6-cff2*adfac1
                ad_cff7=ad_cff7-cff3*adfac1
                ad_cff8=ad_cff8-cff4*adfac1
                ad_dnUdz=ad_dnUdz-                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dnUdx(i,j-1,k1)=ad_dnUdx(i,j-1,k1)+cff1*adfac
                ad_dnUdx(i,j  ,k2)=ad_dnUdx(i,j  ,k2)+cff2*adfac
                ad_dnUdx(i,j-1,k2)=ad_dnUdx(i,j-1,k2)+cff3*adfac
                ad_dnUdx(i,j  ,k1)=ad_dnUdx(i,j  ,k1)+cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j  ,k1)))*        &
!>   &                  tl_dZdx_r(i,j  ,k1)
!>
                ad_dZdx_r(i,j  ,k1)=ad_dZdx_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i,j  ,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j-1,k2)))*        &
!>   &                  tl_dZdx_r(i,j-1,k2)
!>
                ad_dZdx_r(i,j-1,k2)=ad_dZdx_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i,j-1,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j  ,k2)))*        &
!>   &                  tl_dZdx_r(i,j  ,k2)
!>
                ad_dZdx_r(i,j  ,k2)=ad_dZdx_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i,j  ,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j-1,k1)))*        &
!>   &                  tl_dZdx_r(i,j-1,k1)
!>
                ad_dZdx_r(i,j-1,k1)=ad_dZdx_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i,j-1,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
!>   &                  tl_dZde_r(i,j  ,k1)
!>
                ad_dZde_r(i,j  ,k1)=ad_dZde_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j  ,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
!>   &                  tl_dZde_r(i,j-1,k2)
!>
                ad_dZde_r(i,j-1,k2)=ad_dZde_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j-1,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
!>   &                  tl_dZde_r(i,j  ,k2)
!>
                ad_dZde_r(i,j  ,k2)=ad_dZde_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
!>   &                  tl_dZde_r(i,j-1,k1)
!>
                ad_dZde_r(i,j-1,k1)=ad_dZde_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+     &
!>   &                           cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+     &
!>   &                           cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+     &
!>   &                           cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+             &
     &                   cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+             &
     &                   cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+             &
     &                   cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))*            &
     &                  ad_VFse(i,j,k2)
#endif
!>              tl_VFse(i,j,k2)=fac2*                                   &
!>   &                          (tl_cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+  &
!>   &                           tl_cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+  &
!>   &                           tl_cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+  &
!>   &                           tl_cff4*(cff4*dmVdz-dmVde(i,j  ,k1))+  &
!>   &                           cff1*(tl_cff1*dmVdz+cff1*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j-1,k1))+             &
!>   &                           cff2*(tl_cff2*dmVdz+cff2*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j  ,k2))+             &
!>   &                           cff3*(tl_cff3*dmVdz+cff3*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j-1,k2))+             &
!>   &                           cff4*(tl_cff4*dmVdz+cff4*tl_dmVdz-     &
!>   &                                 tl_dmVde(i,j  ,k1)))
!>
                cff=2.0_r8*dmVdz
                adfac=fac2*ad_VFse(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dmVde(i,j-1,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dmVde(i,j  ,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dmVde(i,j-1,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dmVde(i,j  ,k1))*adfac
                ad_dmVdz=ad_dmVdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dmVde(i,j-1,k1)=ad_dmVde(i,j-1,k1)-cff1*adfac
                ad_dmVde(i,j  ,k2)=ad_dmVde(i,j  ,k2)-cff2*adfac
                ad_dmVde(i,j-1,k2)=ad_dmVde(i,j-1,k2)-cff3*adfac
                ad_dmVde(i,j  ,k1)=ad_dmVde(i,j  ,k1)-cff4*adfac
                ad_VFse(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
!>   &                  tl_dZde_r(i,j  ,k1)
!>
                ad_dZde_r(i,j  ,k1)=ad_dZde_r(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j  ,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
!>   &                  tl_dZde_r(i,j-1,k2)
!>
                ad_dZde_r(i,j-1,k2)=ad_dZde_r(i,j-1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i,j-1,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
!>   &                  tl_dZde_r(i,j  ,k2)
!>
                ad_dZde_r(i,j  ,k2)=ad_dZde_r(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
!>   &                  tl_dZde_r(i,j-1,k1)
!>
                ad_dZde_r(i,j-1,k1)=ad_dZde_r(i,j-1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+     &
!>   &                           cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+     &
!>   &                           cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+     &
!>   &                           cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+             &
     &                   cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+             &
     &                   cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+             &
     &                   cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))*            &
     &                  ad_VFsx(i,j,k2)
#endif
!>              tl_VFsx(i,j,k2)=fac1*                                   &
!>   &                          (tl_cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+  &
!>   &                           tl_cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+  &
!>   &                           tl_cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+  &
!>   &                           tl_cff4*(cff4*dnVdz-dnVdx(i+1,j,k1))+  &
!>   &                           cff1*(tl_cff1*dnVdz+cff1*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i  ,j,k1))+             &
!>   &                           cff2*(tl_cff2*dnVdz+cff2*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i+1,j,k2))+             &
!>   &                           cff3*(tl_cff3*dnVdz+cff3*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i  ,j,k2))+             &
!>   &                           cff4*(tl_cff4*dnVdz+cff4*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i+1,j,k1)))
!>
                cff=2.0_r8*dnVdz
                adfac=fac1*ad_VFsx(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dnVdx(i  ,j,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dnVdx(i+1,j,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dnVdx(i  ,j,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dnVdx(i+1,j,k1))*adfac
                ad_dnVdz=ad_dnVdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dnVdx(i  ,j,k1)=ad_dnVdx(i  ,j,k1)-cff1*adfac
                ad_dnVdx(i+1,j,k2)=ad_dnVdx(i+1,j,k2)-cff2*adfac
                ad_dnVdx(i  ,j,k2)=ad_dnVdx(i  ,j,k2)-cff3*adfac
                ad_dnVdx(i+1,j,k1)=ad_dnVdx(i+1,j,k1)-cff4*adfac
                ad_VFsx(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
!>   &                  tl_dZdx_p(i+1,j,k1)
!>
                ad_dZdx_p(i+1,j,k1)=ad_dZdx_p(i+1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
!>   &                  tl_dZdx_p(i  ,j,k2)
!>
                ad_dZdx_p(i  ,j,k2)=ad_dZdx_p(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
!>   &                  tl_dZdx_p(i+1,j,k2)
!>
                ad_dZdx_p(i+1,j,k2)=ad_dZdx_p(i+1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
!>   &                  tl_dZdx_p(i  ,j,k1)
!>
                ad_dZdx_p(i  ,j,k1)=ad_dZdx_p(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
!>              tl_dmVdz=cff*tl_dVdz(i,j,k2)
!>
                ad_dVdz(i,j,k2)=ad_dVdz(i,j,k2)+cff*ad_dmVdz
                ad_dmVdz=0.0_r8
!>              tl_dmUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
!>   &                                tl_dUdz(i+1,j  ,k2)+              &
!>   &                                tl_dUdz(i  ,j-1,k2)+              &
!>   &                                tl_dUdz(i+1,j-1,k2))
!>
                adfac=cff*0.25_r8*ad_dmUdz
                ad_dUdz(i  ,j-1,k2)=ad_dUdz(i  ,j-1,k2)+adfac
                ad_dUdz(i+1,j-1,k2)=ad_dUdz(i+1,j-1,k2)+adfac
                ad_dUdz(i  ,j  ,k2)=ad_dUdz(i  ,j  ,k2)+adfac
                ad_dUdz(i+1,j  ,k2)=ad_dUdz(i+1,j  ,k2)+adfac
                ad_dmUdz=0.0_r8
!
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
!>              tl_dnVdz=cff*tl_dVdz(i,j,k2)
!>
                ad_dVdz(i,j,k2)=ad_dVdz(i,j,k2)+cff*ad_dnVdz
                ad_dnVdz=0.0_r8
!>              tl_dnUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
!>   &                                tl_dUdz(i+1,j  ,k2)+              &
!>   &                                tl_dUdz(i  ,j-1,k2)+              &
!>   &                                tl_dUdz(i+1,j-1,k2))
!>
                adfac=cff*0.25_r8*ad_dnUdz
                ad_dUdz(i  ,j-1,k2)=ad_dUdz(i  ,j-1,k2)+adfac
                ad_dUdz(i+1,j-1,k2)=ad_dUdz(i+1,j-1,k2)+adfac
                ad_dUdz(i  ,j  ,k2)=ad_dUdz(i  ,j  ,k2)+adfac
                ad_dUdz(i+1,j  ,k2)=ad_dUdz(i+1,j  ,k2)+adfac
                ad_dnUdz=0.0_r8
#ifdef VISC_3DCOEF
!>              tl_fac2=tl_cff*om_v(i,j)
!>              tl_fac1=tl_cff*on_v(i,j)
!>
                ad_cff=ad_cff+                                          &
     &                 on_v(i,j)*ad_fac1+om_v(i,j)*ad_fac2
                ad_fac1=0.0_r8
                ad_fac2=0.0_r8
# ifdef UV_U3ADV_SPLIT
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_Vvis3d_r(i,j-1,k  )+tl_Vvis3d_r(i,j,k  )+    &
!>   &                  tl_Vvis3d_r(i,j-1,k+1)+tl_Vvis3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_Vvis3d_r(i,j-1,k  )=ad_Vvis3d_r(i,j-1,k  )+adfac
                ad_Vvis3d_r(i,j  ,k  )=ad_Vvis3d_r(i,j  ,k  )+adfac
                ad_Vvis3d_r(i,j-1,k+1)=ad_Vvis3d_r(i,j-1,k+1)+adfac
                ad_Vvis3d_r(i,j  ,k+1)=ad_Vvis3d_r(i,j  ,k+1)+adfac
                ad_cff=0.0_r8
# else
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_visc3d_r(i,j-1,k  )+tl_visc3d_r(i,j,k  )+    &
!>   &                  tl_visc3d_r(i,j-1,k+1)+tl_visc3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_visc3d_r(i,j-1,k  )=ad_visc3d_r(i,j-1,k  )+adfac
                ad_visc3d_r(i,j  ,k  )=ad_visc3d_r(i,j  ,k  )+adfac
                ad_visc3d_r(i,j-1,k+1)=ad_visc3d_r(i,j-1,k+1)+adfac
                ad_visc3d_r(i,j  ,k+1)=ad_visc3d_r(i,j  ,k+1)+adfac
                ad_cff=0.0_r8
# endif
#endif
              END DO
            END DO
!
            DO j=Jstrm1,Jendp1
              DO i=IstrUm1,Iendp1
#ifdef VISC_3DCOEF
# ifdef UV_U3ADV_SPLIT
                cff=0.125_r8*                                           &
     &              (Uvis3d_r(i-1,j,k  )+Uvis3d_r(i,j,k  )+             &
     &               Uvis3d_r(i-1,j,k+1)+Uvis3d_r(i,j,k+1))
# else
                cff=0.125_r8*                                           &
     &              (visc3d_r(i-1,j,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i-1,j,k+1)+visc3d_r(i,j,k+1))
# endif
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#else
                cff=0.25_r8*(visc4_r(i-1,j)+visc4_r(i,j))
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#endif
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
!
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+     &
!>   &                           cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+     &
!>   &                           cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+     &
!>   &                           cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
!>
                ad_fac2=ad_fac2-                                        &
     &                  (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+             &
     &                   cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+             &
     &                   cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+             &
     &                   cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))*            &
     &                  ad_UFse(i,j,k2)
#endif
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
!>   &                          fac2*                                   &
!>   &                          (tl_cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+  &
!>   &                           tl_cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+  &
!>   &                           tl_cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+  &
!>   &                           tl_cff4*(cff8*dmVdz-dmVde(i  ,j,k1))+  &
!>   &                           cff1*(tl_cff5*dmVdz+cff5*tl_dmVdz-     &
!>   &                                 tl_dmVde(i-1,j,k1))+             &
!>   &                           cff2*(tl_cff6*dmVdz+cff6*tl_dmVdz-     &
!>   &                                 tl_dmVde(i  ,j,k2))+             &
!>   &                           cff3*(tl_cff7*dmVdz+cff7*tl_dmVdz-     &
!>   &                                 tl_dmVde(i-1,j,k2))+             &
!>   &                           cff4*(tl_cff8*dmVdz+cff8*tl_dmVdz-     &
!>   &                                 tl_dmVde(i  ,j,k1)))
!>
                adfac=fac2*ad_UFse(i,j,k2)
                adfac1=adfac*dmVdz
                ad_cff1=ad_cff1-(cff5*dmVdz-dmVde(i-1,j,k1))*adfac
                ad_cff2=ad_cff2-(cff6*dmVdz-dmVde(i  ,j,k2))*adfac
                ad_cff3=ad_cff3-(cff7*dmVdz-dmVde(i-1,j,k2))*adfac
                ad_cff4=ad_cff4-(cff8*dmVdz-dmVde(i  ,j,k1))*adfac
                ad_cff5=ad_cff5-cff1*adfac1
                ad_cff6=ad_cff6-cff2*adfac1
                ad_cff7=ad_cff7-cff3*adfac1
                ad_cff8=ad_cff8-cff4*adfac1
                ad_dmVdz=ad_dmVdz-                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dmVde(i-1,j,k1)=ad_dmVde(i-1,j,k1)+cff1*adfac
                ad_dmVde(i  ,j,k2)=ad_dmVde(i  ,j,k2)+cff2*adfac
                ad_dmVde(i-1,j,k2)=ad_dmVde(i-1,j,k2)+cff3*adfac
                ad_dmVde(i  ,j,k1)=ad_dmVde(i  ,j,k1)+cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_r(i  ,j,k1)))*        &
!>   &                  tl_dZde_r(i  ,j,k1)
!>
                ad_dZde_r(i  ,j,k1)=ad_dZde_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i  ,j,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_r(i-1,j,k2)))*        &
!>   &                  tl_dZde_r(i-1,j,k2)
!>
                ad_dZde_r(i-1,j,k2)=ad_dZde_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_r(i-1,j,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_r(i  ,j,k2)))*        &
!>   &                  tl_dZde_r(i  ,j,k2)
!>
                ad_dZde_r(i  ,j,k2)=ad_dZde_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i  ,j,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_r(i-1,j,k1)))*        &
!>   &                  tl_dZde_r(i-1,j,k1)
!>
                ad_dZde_r(i-1,j,k1)=ad_dZde_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_r(i-1,j,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
!>   &                  tl_dZdx_r(i  ,j,k1)
!>
                ad_dZdx_r(i  ,j,k1)=ad_dZdx_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
!>   &                  tl_dZdx_r(i-1,j,k2)
!>
                ad_dZdx_r(i-1,j,k2)=ad_dZdx_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
!>   &                  tl_dZdx_r(i  ,j,k2)
!>
                ad_dZdx_r(i  ,j,k2)=ad_dZdx_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
!>   &                  tl_dZdx_r(i-1,j,k1)
!>
                ad_dZdx_r(i-1,j,k1)=ad_dZdx_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+     &
!>   &                           cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+     &
!>   &                           cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+     &
!>   &                           cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+             &
     &                   cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+             &
     &                   cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+             &
     &                   cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))*            &
     &                  ad_UFsx(i,j,k2)
#endif
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          fac1*                                   &
!>   &                          (tl_cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+  &
!>   &                           tl_cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+  &
!>   &                           tl_cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+  &
!>   &                           tl_cff4*(cff8*dnVdz-dnVdx(i,j+1,k1))+  &
!>   &                           cff1*(tl_cff5*dnVdz+cff5*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j  ,k1))+             &
!>   &                           cff2*(tl_cff6*dnVdz+cff6*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j+1,k2))+             &
!>   &                           cff3*(tl_cff7*dnVdz+cff7*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j  ,k2))+             &
!>   &                           cff4*(tl_cff8*dnVdz+cff8*tl_dnVdz-     &
!>   &                                 tl_dnVdx(i,j+1,k1)))
!>
                adfac=fac1*ad_UFsx(i,j,k2)
                adfac1=adfac*dnVdz
                ad_cff1=ad_cff1+(cff5*dnVdz-dnVdx(i,j  ,k1))*adfac
                ad_cff2=ad_cff2+(cff6*dnVdz-dnVdx(i,j+1,k2))*adfac
                ad_cff3=ad_cff3+(cff7*dnVdz-dnVdx(i,j  ,k2))*adfac
                ad_cff4=ad_cff4+(cff8*dnVdz-dnVdx(i,j+1,k1))*adfac
                ad_cff5=ad_cff5+cff1*adfac1
                ad_cff6=ad_cff6+cff2*adfac1
                ad_cff7=ad_cff7+cff3*adfac1
                ad_cff8=ad_cff8+cff4*adfac1
                ad_dnVdz=ad_dnVdz+                                      &
     &                   (cff1*cff5+cff2*cff6+cff3*cff7+cff4*cff8)*     &
     &                   adfac
                ad_dnVdx(i,j  ,k1)=ad_dnVdx(i,j  ,k1)-cff1*adfac
                ad_dnVdx(i,j+1,k2)=ad_dnVdx(i,j+1,k2)-cff2*adfac
                ad_dnVdx(i,j  ,k2)=ad_dnVdx(i,j  ,k2)-cff3*adfac
                ad_dnVdx(i,j+1,k1)=ad_dnVdx(i,j+1,k1)-cff4*adfac
!>              tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j+1,k1)))*        &
!>   &                  tl_dZdx_p(i,j+1,k1)
!>
                ad_dZdx_p(i,j+1,k1)=ad_dZdx_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i,j+1,k1)))*   &
     &                              ad_cff8
                ad_cff8=0.0_r8
!>              tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j  ,k2)))*        &
!>   &                  tl_dZdx_p(i,j  ,k2)
!>
                ad_dZdx_p(i,j  ,k2)=ad_dZdx_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_p(i,j  ,k2)))*   &
     &                              ad_cff7
                ad_cff7=0.0_r8
!>              tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j+1,k2)))*        &
!>   &                  tl_dZdx_p(i,j+1,k2)
!>
                ad_dZdx_p(i,j+1,k2)=ad_dZdx_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i,j+1,k2)))*   &
     &                              ad_cff6
                ad_cff6=0.0_r8
!>              tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j  ,k1)))*        &
!>   &                  tl_dZdx_p(i,j  ,k1)
!>
                ad_dZdx_p(i,j  ,k1)=ad_dZdx_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_p(i,j  ,k1)))*   &
     &                              ad_cff5
                ad_cff5=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
!>   &                  tl_dZde_p(i,j+1,k1)
!>
                ad_dZde_p(i,j+1,k1)=ad_dZde_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j+1,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
!>   &                  tl_dZde_p(i,j  ,k2)
!>
                ad_dZde_p(i,j  ,k2)=ad_dZde_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j  ,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
!>   &                  tl_dZde_p(i,j+1,k2)
!>
                ad_dZde_p(i,j+1,k2)=ad_dZde_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
!>   &                  tl_dZde_p(i,j  ,k1)
!>
                ad_dZde_p(i,j  ,k1)=ad_dZde_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFse(i,j,k2)=tl_UFse(i,j,k2)+
!>   &                          tl_fac2*                                &
!>   &                          (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+     &
!>   &                           cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+     &
!>   &                           cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+     &
!>   &                           cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
!>
                ad_fac2=ad_fac2+                                        &
     &                  (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+             &
     &                   cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+             &
     &                   cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+             &
     &                   cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))*            &
     &                  ad_UFse(i,j,k2)
#endif
!>              tl_UFse(i,j,k2)=fac2*                                   &
!>   &                          (tl_cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+  &
!>   &                           tl_cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+  &
!>   &                           tl_cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+  &
!>   &                           tl_cff4*(cff4*dmUdz-dmUde(i,j+1,k1))+  &
!>   &                           cff1*(tl_cff1*dmUdz+cff1*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j  ,k1))+             &
!>   &                           cff2*(tl_cff2*dmUdz+cff2*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j+1,k2))+             &
!>   &                           cff3*(tl_cff3*dmUdz+cff3*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j  ,k2))+             &
!>   &                           cff4*(tl_cff4*dmUdz+cff4*tl_dmUdz-     &
!>   &                                 tl_dmUde(i,j+1,k1)))
!>
                cff=2.0_r8*dmUdz
                adfac=fac2*ad_UFse(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dmUde(i,j  ,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dmUde(i,j+1,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dmUde(i,j  ,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dmUde(i,j+1,k1))*adfac
                ad_dmUdz=ad_dmUdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dmUde(i,j  ,k1)=ad_dmUde(i,j  ,k1)-cff1*adfac
                ad_dmUde(i,j+1,k2)=ad_dmUde(i,j+1,k2)-cff2*adfac
                ad_dmUde(i,j  ,k2)=ad_dmUde(i,j  ,k2)-cff3*adfac
                ad_dmUde(i,j+1,k1)=ad_dmUde(i,j+1,k1)-cff4*adfac
                ad_UFse(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
!>   &                  tl_dZde_p(i,j+1,k1)
!>
                ad_dZde_p(i,j+1,k1)=ad_dZde_p(i,j+1,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j+1,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
!>   &                  tl_dZde_p(i,j  ,k2)
!>
                ad_dZde_p(i,j  ,k2)=ad_dZde_p(i,j  ,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZde_p(i,j  ,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
!>   &                  tl_dZde_p(i,j+1,k2)
!>
                ad_dZde_p(i,j+1,k2)=ad_dZde_p(i,j+1,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
!>   &                  tl_dZde_p(i,j  ,k1)
!>
                ad_dZde_p(i,j  ,k1)=ad_dZde_p(i,j  ,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
!>              tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
!>   &                          tl_fac1*                                &
!>   &                          (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+     &
!>   &                           cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+     &
!>   &                           cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+     &
!>   &                           cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
!>
                ad_fac1=ad_fac1+                                        &
     &                  (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+             &
     &                   cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+             &
     &                   cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+             &
     &                   cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))*            &
     &                  ad_UFsx(i,j,k2)
#endif
!>              tl_UFsx(i,j,k2)=fac1*                                   &
!>   &                          (tl_cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+  &
!>   &                           tl_cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+  &
!>   &                           tl_cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+  &
!>   &                           tl_cff4*(cff4*dnUdz-dnUdx(i  ,j,k1))+  &
!>   &                           cff1*(tl_cff1*dnUdz+cff1*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i-1,j,k1))+             &
!>   &                           cff2*(tl_cff2*dnUdz+cff2*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i  ,j,k2))+             &
!>   &                           cff3*(tl_cff3*dnUdz+cff3*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i-1,j,k2))+             &
!>   &                           cff4*(tl_cff4*dnUdz+cff4*tl_dnUdz-     &
!>   &                                 tl_dnUdx(i  ,j,k1)))
!>
                cff=2.0_r8*dnUdz
                adfac=fac1*ad_UFsx(i,j,k2)
                ad_cff1=ad_cff1+(cff1*cff-dnUdx(i-1,j,k1))*adfac
                ad_cff2=ad_cff2+(cff2*cff-dnUdx(i  ,j,k2))*adfac
                ad_cff3=ad_cff3+(cff3*cff-dnUdx(i-1,j,k2))*adfac
                ad_cff4=ad_cff4+(cff4*cff-dnUdx(i  ,j,k1))*adfac
                ad_dnUdz=ad_dnUdz+                                      &
     &                   (cff1*cff1+cff2*cff2+cff3*cff3+cff4*cff4)*     &
     &                   adfac
                ad_dnUdx(i-1,j,k1)=ad_dnUdx(i-1,j,k1)-cff1*adfac
                ad_dnUdx(i  ,j,k2)=ad_dnUdx(i  ,j,k2)-cff2*adfac
                ad_dnUdx(i-1,j,k2)=ad_dnUdx(i-1,j,k2)-cff3*adfac
                ad_dnUdx(i  ,j,k1)=ad_dnUdx(i  ,j,k1)-cff4*adfac
                ad_UFsx(i,j,k2)=0.0_r8
!>              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
!>   &                  tl_dZdx_r(i  ,j,k1)
!>
                ad_dZdx_r(i  ,j,k1)=ad_dZdx_r(i  ,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*   &
     &                              ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
!>   &                  tl_dZdx_r(i-1,j,k2)
!>
                ad_dZdx_r(i-1,j,k2)=ad_dZdx_r(i-1,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*   &
     &                              ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
!>   &                  tl_dZdx_r(i  ,j,k2)
!>
                ad_dZdx_r(i  ,j,k2)=ad_dZdx_r(i  ,j,k2)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*   &
     &                              ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
!>   &                  tl_dZdx_r(i-1,j,k1)
!>
                ad_dZdx_r(i-1,j,k1)=ad_dZdx_r(i-1,j,k1)+                &
     &                              (0.5_r8+                            &
     &                               SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*   &
     &                              ad_cff1
                ad_cff1=0.0_r8
!
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
!>              tl_dmVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
!>   &                                tl_dVdz(i  ,j+1,k2)+              &
!>   &                                tl_dVdz(i-1,j  ,k2)+              &
!>   &                                tl_dVdz(i  ,j  ,k2))
!>
                adfac=cff*0.25_r8*ad_dmVdz
                ad_dVdz(i-1,j  ,k2)=ad_dVdz(i-1,j  ,k2)+adfac
                ad_dVdz(i  ,j  ,k2)=ad_dVdz(i  ,j  ,k2)+adfac
                ad_dVdz(i-1,j+1,k2)=ad_dVdz(i-1,j+1,k2)+adfac
                ad_dVdz(i  ,j+1,k2)=ad_dVdz(i  ,j+1,k2)+adfac
                ad_dmVdz=0.0_r8
!>              tl_dmUdz=cff*tl_dUdz(i,j,k2)
!>
                ad_dUdz(i,j,k2)=ad_dUdz(i,j,k2)+cff*ad_dmUdz
                ad_dmUdz=0.0_r8
!
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
!>              tl_dnVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
!>   &                                tl_dVdz(i  ,j+1,k2)+              &
!>   &                                tl_dVdz(i-1,j  ,k2)+              &
!>   &                                tl_dVdz(i  ,j  ,k2))
!>
                adfac=cff*0.25_r8*ad_dnVdz
                ad_dVdz(i-1,j  ,k2)=ad_dVdz(i-1,j  ,k2)+adfac
                ad_dVdz(i  ,j  ,k2)=ad_dVdz(i  ,j  ,k2)+adfac
                ad_dVdz(i-1,j+1,k2)=ad_dVdz(i-1,j+1,k2)+adfac
                ad_dVdz(i  ,j+1,k2)=ad_dVdz(i  ,j+1,k2)+adfac
                ad_dnVdz=0.0_r8
!>              tl_dnUdz=cff*tl_dUdz(i,j,k2)
!>
                ad_dUdz(i,j,k2)=ad_dUdz(i,j,k2)+cff*ad_dnUdz
                ad_dnUdz=0.0_r8
#ifdef VISC_3DCOEF
!>              tl_fac2=tl_cff*om_u(i,j)
!>              tl_fac1=tl_cff*on_u(i,j)
!>
                ad_cff=ad_cff+                                          &
     &                 on_u(i,j)*ad_fac1+om_u(i,j)*ad_fac2
                ad_fac1=0.0_r8
                ad_fac2=0.0_r8
# ifdef UV_U3ADV_SPLIT
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_Uvis3d_r(i-1,j,k  )+tl_Uvis3d_r(i,j,k  )+    &
!>   &                  tl_Uvis3d_r(i-1,j,k+1)+tl_Uvis3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_Uvis3d_r(i-1,j,k  )=ad_Uvis3d_r(i-1,j,k  )+adfac
                ad_Uvis3d_r(i  ,j,k  )=ad_Uvis3d_r(i  ,j,k  )+adfac
                ad_Uvis3d_r(i-1,j,k+1)=ad_Uvis3d_r(i-1,j,k+1)+adfac
                ad_Uvis3d_r(i  ,j,k+1)=ad_Uvis3d_r(i  ,j,k+1)+adfac
                ad_cff=0.0_r8
# else
!>              tl_cff=0.125_r8*                                        &
!>   &                 (tl_visc3d_r(i-1,j,k  )+tl_visc3d_r(i,j,k  )+    &
!>   &                  tl_visc3d_r(i-1,j,k+1)+tl_visc3d_r(i,j,k+1))
!>
                adfac=0.125_r8*ad_cff
                ad_visc3d_r(i-1,j,k  )=ad_visc3d_r(i-1,j,k  )+adfac
                ad_visc3d_r(i  ,j,k  )=ad_visc3d_r(i  ,j,k  )+adfac
                ad_visc3d_r(i-1,j,k+1)=ad_visc3d_r(i-1,j,k+1)+adfac
                ad_visc3d_r(i  ,j,k+1)=ad_visc3d_r(i  ,j,k+1)+adfac
                ad_cff=0.0_r8
# endif
#endif
              END DO
            END DO
          END IF
        END IF
!
!  Compute adjoint components of the rotated viscous flux (m^4 s-^3/2)
!  along geopotential surfaces in the XI- and ETA-directions.
!
        IF (k.gt.0) THEN
          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              pm_p=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                      &
     &                      pm(i  ,j-1)+pm(i  ,j))
              pn_p=0.25_r8*(pn(i-1,j-1)+pn(i-1,j)+                      &
     &                      pn(i  ,j-1)+pn(i  ,j))
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
              cff3=MIN(dZde_p(i,j,k1),0.0_r8)
              cff4=MAX(dZde_p(i,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
              cff=on_p(i,j)*(dnVdx(i,j,k1)-                             &
     &                       0.5_r8*pn_p*                               &
     &                       (cff1*(dVdz(i-1,j,k1)+                     &
     &                              dVdz(i  ,j,k2))+                    &
     &                        cff2*(dVdz(i-1,j,k2)+                     &
     &                              dVdz(i  ,j,k1))))+                  &
     &            om_p(i,j)*(dmUde(i,j,k1)-                             &
     &                       0.5_r8*pm_p*                               &
     &                       (cff3*(dUdz(i,j-1,k1)+                     &
     &                              dUdz(i,j  ,k2))+                    &
     &                        cff4*(dUdz(i,j-1,k2)+                     &
     &                              dUdz(i,j  ,k1))))
# ifdef MASKING
              cff=cff*pmask(i,j)
# endif
# ifdef UV_U3ADV_SPLIT
              Uvis_p=0.25_r8*                                           &
     &               (Uvis3d_r(i-1,j-1,k)+Uvis3d_r(i-1,j,k)+            &
     &                Uvis3d_r(i  ,j-1,k)+Uvis3d_r(i  ,j,k))
              Vvis_p=0.25_r8*                                           &
     &               (Vvis3d_r(i-1,j-1,k)+Vvis3d_r(i-1,j,k)+            &
     &                Vvis3d_r(i  ,j-1,k)+Vvis3d_r(i  ,j,k))
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*                          &
!>   &                    (tl_Vvis_p*cff+Vvis_p*tl_cff)
!>
              adfac=on_p(i,j)*on_p(i,j)*ad_VFx(i,j)
              ad_cff=ad_cff+Vvis_p*adfac
              ad_Vvis_p=ad_Vvis_p+cff*adfac
              ad_VFx(i,j)=0.0_r8
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*                          &
!>   &                    (tl_Uvis_p*cff+Uvis_p*tl_cff)
!>
              adfac=om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
              ad_cff=ad_cff+Uvis_p*adfac
              ad_Uvis_p=ad_Uvis_p+cff*adfac
              ad_UFe(i,j)=0.0_r8
!>            tl_Vvis_p=0.25_r8*                                        &
!>   &                  (tl_Vvis3d_r(i-1,j-1,k)+tl_Vvis3d_r(i-1,j,k)+   &
!>   &                   tl_Vvis3d_r(i  ,j-1,k)+tl_Vvis3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_Vvis_p
              ad_Vvis3d_r(i-1,j-1,k)=ad_Vvis3d_r(i-1,j-1,k)+adfac
              ad_Vvis3d_r(i-1,j  ,k)=ad_Vvis3d_r(i-1,j  ,k)+adfac
              ad_Vvis3d_r(i  ,j-1,k)=ad_Vvis3d_r(i  ,j-1,k)+adfac
              ad_Vvis3d_r(i  ,j  ,k)=ad_Vvis3d_r(i  ,j  ,k)+adfac
              ad_Vvis_p=0.0_r8
!>            tl_Uvis_p=0.25_r8*                                        &
!>   &                  (tl_Uvis3d_r(i-1,j-1,k)+tl_Uvis3d_r(i-1,j,k)+   &
!>   &                   tl_Uvis3d_r(i  ,j-1,k)+tl_Uvis3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_Uvis_p
              ad_Uvis3d_r(i-1,j-1,k)=ad_Uvis3d_r(i-1,j-1,k)+adfac
              ad_Uvis3d_r(i-1,j  ,k)=ad_Uvis3d_r(i-1,j  ,k)+adfac
              ad_Uvis3d_r(i  ,j-1,k)=ad_Uvis3d_r(i  ,j-1,k)+adfac
              ad_Uvis3d_r(i  ,j  ,k)=ad_Uvis3d_r(i  ,j  ,k)+adfac
              ad_Uvis_p=0.0_r8
# else
              visc_p=0.25_r8*                                           &
     &               (visc3d_r(i-1,j-1,k)+visc3d_r(i-1,j,k)+            &
     &                visc3d_r(i  ,j-1,k)+visc3d_r(i  ,j,k))
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*                          &
!>   &                    (tl_visc_p*cff+visc_p*tl_cff)
!>
              adfac=on_p(i,j)*on_p(i,j)*ad_VFx(i,j)
              ad_cff=ad_cff+visc_p*adfac
              ad_visc_p=ad_visc_p+cff*adfac
              ad_VFx(i,j)=0.0_r8
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*                          &
!>   &                    (tl_visc_p*cff+visc_p*tl_cff)
!>
              adfac=om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
              ad_cff=ad_cff+visc_p*adfac
              ad_visc_p=ad_visc_p+cff*adfac
              ad_UFe(i,j)=0.0_r8
!>            tl_visc_p=0.25_r8*                                        &
!>   &                  (tl_visc3d_r(i-1,j-1,k)+tl_visc3d_r(i-1,j,k)+   &
!>   &                   tl_visc3d_r(i  ,j-1,k)+tl_visc3d_r(i  ,j,k))
!>
              adfac=0.25_r8*ad_visc_p
              ad_visc3d_r(i-1,j-1,k)=ad_visc3d_r(i-1,j-1,k)+adfac
              ad_visc3d_r(i  ,j-1,k)=ad_visc3d_r(i  ,j-1,k)+adfac
              ad_visc3d_r(i-1,j  ,k)=ad_visc3d_r(i-1,j  ,k)+adfac
              ad_visc3d_r(i  ,j  ,k)=ad_visc3d_r(i  ,j  ,k)+adfac
              ad_visc_p=0.0_r8
# endif
#else
!>            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*visc4_p(i,j)*tl_cff
!>            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*visc4_p(i,j)*tl_cff
!>
              ad_cff=ad_cff+                                            &
     &               on_p(i,j)*on_p(i,j)*visc4_p(i,j)*ad_VFx(i,j)+      &
     &               om_p(i,j)*om_p(i,j)*visc4_p(i,j)*ad_UFe(i,j)
              ad_VFx(i,j)=0.0_r8
              ad_UFe(i,j)=0.0_r8
#endif
#ifdef MASKING
!>            tl_cff=tl_cff*pmask(i,j)
!>
              ad_cff=ad_cff*pmask(i,j)
#endif
!>            tl_cff=on_p(i,j)*(tl_dnVdx(i,j,k1)-                       &
!>   &                          0.5_r8*pn_p*                            &
!>   &                          (tl_cff1*(dVdz(i-1,j,k1)+               &
!>   &                                    dVdz(i  ,j,k2))+              &
!>   &                           cff1*(tl_dVdz(i-1,j,k1)+               &
!>   &                                 tl_dVdz(i  ,j,k2))+              &
!>   &                           tl_cff2*(dVdz(i-1,j,k2)+               &
!>   &                                    dVdz(i  ,j,k1))+              &
!>   &                           cff2*(tl_dVdz(i-1,j,k2)+               &
!>   &                                 tl_dVdz(i  ,j,k1))))+            &
!>   &               om_p(i,j)*(tl_dmUde(i,j,k1)-                       &
!>   &                          0.5_r8*pm_p*                            &
!>   &                          (tl_cff3*(dUdz(i,j-1,k1)+               &
!>   &                                    dUdz(i,j  ,k2))+              &
!>   &                           cff3*(tl_dUdz(i,j-1,k1)+               &
!>   &                                 tl_dUdz(i,j  ,k2))+              &
!>   &                           tl_cff4*(dUdz(i,j-1,k2)+               &
!>   &                                    dUdz(i,j  ,k1))+              &
!>   &                           cff4*(tl_dUdz(i,j-1,k2)+               &
!>   &                                 tl_dUdz(i,j  ,k1))))
!>
              adfac1=on_p(i,j)*ad_cff
              adfac2=adfac1*0.5_r8*pn_p
              adfac3=om_p(i,j)*ad_cff
              adfac4=adfac3*0.5_r8*pm_p
              ad_dnVdx(i,j,k1)=ad_dnVdx(i,j,k1)+adfac1
              ad_cff1=ad_cff1-                                          &
     &                (dVdz(i-1,j,k1)+dVdz(i  ,j,k2))*adfac2
              ad_cff2=ad_cff2-                                          &
     &                (dVdz(i-1,j,k2)+dVdz(i  ,j,k1))*adfac2
              ad_dVdz(i-1,j,k1)=ad_dVdz(i-1,j,k1)-cff1*adfac2
              ad_dVdz(i-1,j,k2)=ad_dVdz(i-1,j,k2)-cff2*adfac2
              ad_dVdz(i  ,j,k1)=ad_dVdz(i  ,j,k1)-cff2*adfac2
              ad_dVdz(i  ,j,k2)=ad_dVdz(i  ,j,k2)-cff1*adfac2
              ad_dmUde(i,j,k1)=ad_dmUde(i,j,k1)+adfac3
              ad_cff3=ad_cff3-                                          &
     &                (dUdz(i,j-1,k1)+dUdz(i,j  ,k2))*adfac4
              ad_cff4=ad_cff4-                                          &
     &                (dUdz(i,j-1,k2)+dUdz(i,j  ,k1))*adfac4
              ad_dUdz(i,j-1,k1)=ad_dUdz(i,j-1,k1)-cff3*adfac4
              ad_dUdz(i,j-1,k2)=ad_dUdz(i,j-1,k2)-cff4*adfac4
              ad_dUdz(i,j  ,k1)=ad_dUdz(i,j  ,k1)-cff4*adfac4
              ad_dUdz(i,j  ,k2)=ad_dUdz(i,j  ,k2)-cff3*adfac4
              ad_cff=0.0_r8
!>            tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j,k1)))*            &
!>   &                tl_dZde_p(i,j,k1)
!>            tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j,k1)))*            &
!>   &                tl_dZde_p(i,j,k1)
!>
              ad_dZde_p(i,j,k1)=ad_dZde_p(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZde_p(i,j,k1)))*         &
     &                          ad_cff4+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZde_p(i,j,k1)))*         &
     &                          ad_cff3
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j,k1)))*            &
!>   &                tl_dZdx_p(i,j,k1)
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j,k1)))*            &
!>   &                tl_dZdx_p(i,j,k1)
!>
              ad_dZdx_p(i,j,k1)=ad_dZdx_p(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZdx_p(i,j,k1)))*         &
     &                          ad_cff2+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZdx_p(i,j,k1)))*         &
     &                          ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO
!
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
              cff3=MIN(dZde_r(i,j,k1),0.0_r8)
              cff4=MAX(dZde_r(i,j,k1),0.0_r8)
#ifdef VISC_3DCOEF
              cff=on_r(i,j)*(dnUdx(i,j,k1)-                             &
     &                       0.5_r8*pn(i,j)*                            &
     &                       (cff1*(dUdz(i  ,j,k1)+                     &
     &                              dUdz(i+1,j,k2))+                    &
     &                        cff2*(dUdz(i  ,j,k2)+                     &
     &                              dUdz(i+1,j,k1))))-                  &
     &            om_r(i,j)*(dmVde(i,j,k1)-                             &
     &                       0.5_r8*pm(i,j)*                            &
     &                       (cff3*(dVdz(i,j  ,k1)+                     &
     &                              dVdz(i,j+1,k2))+                    &
     &                        cff4*(dVdz(i,j  ,k2)+                     &
     &                              dVdz(i,j+1,k1))))
# ifdef MASKING
              cff=cff*rmask(i,j)
# endif
# ifdef UV_U3ADV_SPLIT
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*                          &
!>   &                    (tl_Vvis3d_r(i,j,k)*cff+                      &
!>   &                     Vvis3d_r(i,j,k)*tl_cff)
!>
              adfac=om_r(i,j)*om_r(i,j)*ad_VFe(i,j)
              ad_cff=ad_cff+Vvis3d_r(i,j,k)*adfac
              ad_Vvis3d_r(i,j,k)=ad_Vvis3d_r(i,j,k)+cff*adfac
              ad_VFe(i,j)=0.0_r8
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*
!>   &                    (tl_Uvis3d_r(i,j,k)*cff+                      &
!>   &                     Uvis3d_r(i,j,k)*tl_cff)
!>
              adfac=on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
              ad_cff=ad_cff+Uvis3d_r(i,j,k)*adfac
              ad_Uvis3d_r(i,j,k)=ad_Uvis3d_r(i,j,k)+cff*adfac
              ad_UFx(i,j)=0.0_r8
# else
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*                          &
!>   &                    (tl_visc3d_r(i,j,k)*cff+                      &
!>   &                     visc3d_r(i,j,k)*tl_cff)
!>
              adfac=om_r(i,j)*om_r(i,j)*ad_VFe(i,j)
              ad_cff=ad_cff+visc3d_r(i,j,k)*adfac
              ad_visc3d_r(i,j,k)=ad_visc3d_r(i,j,k)+cff*adfac
              ad_VFe(i,j)=0.0_r8
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*                          &
!>   &                    (tl_visc3d_r(i,j,k)*cff+                      &
!>   &                     visc3d_r(i,j,k)*tl_cff)
!>
              adfac=on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
              ad_cff=ad_cff+visc3d_r(i,j,k)*adfac
              ad_visc3d_r(i,j,k)=ad_visc3d_r(i,j,k)+cff*adfac
              ad_UFx(i,j)=0.0_r8
# endif
#else
!>            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*visc4_r(i,j)*tl_cff
!>            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*visc4_r(i,j)*tl_cff
!>
              ad_cff=ad_cff+                                            &
     &               om_r(i,j)*om_r(i,j)*visc4_r(i,j)*ad_VFe(i,j)+      &
     &               on_r(i,j)*on_r(i,j)*visc4_r(i,j)*ad_UFx(i,j)
              ad_VFe(i,j)=0.0_r8
              ad_UFx(i,j)=0.0_r8
#endif
#ifdef MASKING
!>            tl_cff=tl_cff*rmask(i,j)
!>
              ad_cff=ad_cff*rmask(i,j)
#endif
!>            tl_cff=on_r(i,j)*(tl_dnUdx(i,j,k1)-                       &
!>   &                          0.5_r8*pn(i,j)*                         &
!>   &                          (tl_cff1*(dUdz(i  ,j,k1)+               &
!>   &                                    dUdz(i+1,j,k2))+              &
!>   &                           cff1*(tl_dUdz(i  ,j,k1)+               &
!>   &                                 tl_dUdz(i+1,j,k2))+              &
!>   &                           tl_cff2*(dUdz(i  ,j,k2)+               &
!>   &                                    dUdz(i+1,j,k1))+              &
!>   &                           cff2*(tl_dUdz(i  ,j,k2)+               &
!>   &                                 tl_dUdz(i+1,j,k1))))-            &
!>   &               om_r(i,j)*(tl_dmVde(i,j,k1)-                       &
!>   &                          0.5_r8*pm(i,j)*                         &
!>   &                          (tl_cff3*(dVdz(i,j  ,k1)+               &
!>   &                                    dVdz(i,j+1,k2))+              &
!>   &                           cff3*(tl_dVdz(i,j  ,k1)+               &
!>   &                                 tl_dVdz(i,j+1,k2))+              &
!>   &                           tl_cff4*(dVdz(i,j  ,k2)+               &
!>   &                                    dVdz(i,j+1,k1))+              &
!>   &                           cff4*(tl_dVdz(i,j  ,k2)+               &
!>   &                                 tl_dVdz(i,j+1,k1))))
!>
              adfac1=on_r(i,j)*ad_cff
              adfac2=adfac1*0.5_r8*pn(i,j)
              adfac3=om_r(i,j)*ad_cff
              adfac4=adfac3*0.5_r8*pm(i,j)
              ad_dnUdx(i,j,k1)=ad_dnUdx(i,j,k1)+adfac1
              ad_cff1=ad_cff1-                                          &
     &                (dUdz(i  ,j,k1)+dUdz(i+1,j,k2))*adfac2
              ad_cff2=ad_cff2-                                          &
     &                (dUdz(i  ,j,k2)+dUdz(i+1,j,k1))*adfac2
              ad_dUdz(i  ,j,k1)=ad_dUdz(i  ,j,k1)-cff1*adfac2
              ad_dUdz(i  ,j,k2)=ad_dUdz(i  ,j,k2)-cff2*adfac2
              ad_dUdz(i+1,j,k1)=ad_dUdz(i+1,j,k1)-cff2*adfac2
              ad_dUdz(i+1,j,k2)=ad_dUdz(i+1,j,k2)-cff1*adfac2
              ad_dmVde(i,j,k1)=ad_dmVde(i,j,k1)-adfac3
              ad_cff3=ad_cff3+                                          &
     &                (dVdz(i,j  ,k1)+dVdz(i,j+1,k2))*adfac4
              ad_cff4=ad_cff4+                                          &
     &                (dVdz(i,j  ,k2)+dVdz(i,j+1,k1))*adfac4
              ad_dVdz(i,j  ,k1)=ad_dVdz(i,j  ,k1)+cff3*adfac4
              ad_dVdz(i,j  ,k2)=ad_dVdz(i,j  ,k2)+cff4*adfac4
              ad_dVdz(i,j+1,k1)=ad_dVdz(i,j+1,k1)+cff4*adfac4
              ad_dVdz(i,j+1,k2)=ad_dVdz(i,j+1,k2)+cff3*adfac4
              ad_cff=0.0_r8
!>            tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j,k1)))*            &
!>   &                tl_dZde_r(i,j,k1)
!>            tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j,k1)))*            &
!>   &                tl_dZde_r(i,j,k1)
!>
              ad_dZde_r(i,j,k1)=ad_dZde_r(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZde_r(i,j,k1)))*         &
     &                          ad_cff4+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZde_r(i,j,k1)))*         &
     &                          ad_cff3
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j,k1)))*            &
!>   &                tl_dZdx_r(i,j,k1)
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j,k1)))*            &
!>   &                tl_dZdx_r(i,j,k1)
!>
              ad_dZdx_r(i,j,k1)=ad_dZdx_r(i,j,k1)+                      &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8, dZdx_r(i,j,k1)))*         &
     &                          ad_cff2+                                &
     &                          (0.5_r8+                                &
     &                           SIGN(0.5_r8,-dZdx_r(i,j,k1)))*         &
     &                          ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO
        END IF
!
!  Compute momentum horizontal (1/m/s) and vertical (1/s) adjoint
!  gradients.
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=JstrVm1,Jendp1
            DO i=Istrm1,Iendp1
!>            tl_VFse(i,j,k2)=0.0_r8
!>
              ad_VFse(i,j,k2)=0.0_r8
!>            tl_VFsx(i,j,k2)=0.0_r8
!>
              ad_VFsx(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstrm1,Jendp1
            DO i=IstrUm1,Iendp1
!>            tl_UFse(i,j,k2)=0.0_r8
!>
              ad_UFse(i,j,k2)=0.0_r8
!>            tl_UFsx(i,j,k2)=0.0_r8
!>
              ad_UFsx(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
!>            tl_dVdz(i,j,k2)=0.0_r8
!>
              ad_dVdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
!>            tl_dUdz(i,j,k2)=0.0_r8
!>
              ad_dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-                       &
     &                            z_r(i,j-1,k  )+                       &
     &                            z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j  ,k  )))
!>            tl_dVdz(i,j,k2)=tl_cff*(v(i,j,k+1,nrhs)-                  &
!>   &                                v(i,j,k  ,nrhs))+                 &
!>   &                        cff*(tl_v(i,j,k+1,nrhs)-                  &
!>   &                             tl_v(i,j,k  ,nrhs))
!>
              adfac=cff*ad_dVdz(i,j,k2)
              ad_v(i,j,k  ,nrhs)=ad_v(i,j,k  ,nrhs)-adfac
              ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+adfac
              ad_cff=ad_cff+(v(i,j,k+1,nrhs)-                           &
     &                       v(i,j,k  ,nrhs))*ad_dVdz(i,j,k2)
              ad_dVdz(i,j,k2)=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(tl_z_r(i,j-1,k+1)-               &
!>   &                                 tl_z_r(i,j-1,k  )+               &
!>   &                                 tl_z_r(i,j  ,k+1)-               &
!>   &                                 tl_z_r(i,j  ,k  )))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+adfac
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO

          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-                       &
     &                            z_r(i-1,j,k  )+                       &
     &                            z_r(i  ,j,k+1)-                       &
     &                            z_r(i  ,j,k  )))
!>            tl_dUdz(i,j,k2)=tl_cff*(u(i,j,k+1,nrhs)-                  &
!>   &                                u(i,j,k  ,nrhs))+                 &
!>   &                        cff*(tl_u(i,j,k+1,nrhs)-                  &
!>   &                             tl_u(i,j,k  ,nrhs))
!>
              adfac=cff*ad_dUdz(i,j,k2)
              ad_u(i,j,k  ,nrhs)=ad_u(i,j,k  ,nrhs)-adfac
              ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+adfac
              ad_cff=ad_cff+(u(i,j,k+1,nrhs)-                           &
     &                       u(i,j,k  ,nrhs))*ad_dUdz(i,j,k2)
              ad_dUdz(i,j,k2)=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*((tl_z_r(i-1,j,k+1)-              &
!>   &                                  tl_z_r(i-1,j,k  )+              &
!>   &                                  tl_z_r(i  ,j,k+1)-              &
!>   &                                  tl_z_r(i  ,j,k  )))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac
              ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+adfac
              ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-adfac
              ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO
        END IF

        IF (k.lt.N(ng)) THEN
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff=0.5_r8*pn(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
!>            tl_dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*              &
!>   &                              tl_v(i,j+1,k+1,nrhs)-               &
!>   &                              (pm(i,j-1)+pm(i,j  ))*              &
!>   &                              tl_v(i,j  ,k+1,nrhs))
!>
              adfac=cff*ad_dmVde(i,j,k2)
              ad_v(i,j  ,k+1,nrhs)=ad_v(i,j  ,k+1,nrhs)-                &
     &                             (pm(i,j-1)+pm(i,j  ))*adfac
              ad_v(i,j+1,k+1,nrhs)=ad_v(i,j+1,k+1,nrhs)+                &
     &                             (pm(i,j  )+pm(i,j+1))*adfac
              ad_dmVde(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                      &
     &                      pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
!>            tl_dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*            &
!>   &                              tl_v(i  ,j,k+1,nrhs)-               &
!>   &                              (pn(i-1,j-1)+pn(i-1,j))*            &
!>   &                              tl_v(i-1,j,k+1,nrhs))
!>
              adfac=cff*ad_dnVdx(i,j,k2)
              ad_v(i-1,j,k+1,nrhs)=ad_v(i-1,j,k+1,nrhs)-                &
     &                             (pn(i-1,j-1)+pn(i-1,j))*adfac
              ad_v(i  ,j,k+1,nrhs)=ad_v(i  ,j,k+1,nrhs)+                &
     &                             (pn(i  ,j-1)+pn(i  ,j))*adfac
              ad_dnVdx(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
              cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                      &
     &                      pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
!>            tl_dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*            &
!>   &                              tl_u(i,j  ,k+1,nrhs)-               &
!>   &                              (pm(i-1,j-1)+pm(i,j-1))*            &
!>   &                              tl_u(i,j-1,k+1,nrhs))
!>
              adfac=cff*ad_dmUde(i,j,k2)
              ad_u(i,j-1,k+1,nrhs)=ad_u(i,j-1,k+1,nrhs)-                &
     &                             (pm(i-1,j-1)+pm(i,j-1))*adfac
              ad_u(i,j  ,k+1,nrhs)=ad_u(i,j  ,k+1,nrhs)+                &
     &                             (pm(i-1,j  )+pm(i,j  ))*adfac
              ad_dmUde(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
              cff=0.5_r8*pm(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
!>            tl_dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*              &
!>   &                              tl_u(i+1,j,k+1,nrhs)-               &
!>   &                              (pn(i-1,j)+pn(i  ,j))*              &
!>   &                              tl_u(i  ,j,k+1,nrhs))
!>
              adfac=cff*ad_dnUdx(i,j,k2)
              ad_u(i  ,j,k+1,nrhs)=ad_u(i  ,j,k+1,nrhs)-                &
     &                             (pn(i-1,j)+pn(i  ,j))*adfac
              ad_u(i+1,j,k+1,nrhs)=ad_u(i+1,j,k+1,nrhs)+                &
     &                             (pn(i  ,j)+pn(i+1,j))*adfac
              ad_dnUdx(i,j,k2)=0.0_r8
            END DO
          END DO
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=JstrVm2,Jendp1
            DO i=IstrUm2,Iendp1
!>            tl_dZde_r(i,j,k2)=0.5_r8*(tl_VFe(i,j  )+                  &
!>   &                                  tl_VFe(i,j+1))
!>
              adfac=0.5_r8*ad_dZde_r(i,j,k2)
              ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
              ad_VFe(i,j+1)=ad_VFe(i,j+1)+adfac
              ad_dZde_r(i,j,k2)=0.0_r8
!>            tl_dZdx_r(i,j,k2)=0.5_r8*(tl_UFx(i  ,j)+                  &
!>   &                                  tl_UFx(i+1,j))
!>
              adfac=0.5_r8*ad_dZdx_r(i,j,k2)
              ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
              ad_UFx(i+1,j)=ad_UFx(i+1,j)+adfac
              ad_dZdx_r(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstrm1,Jendp2
            DO i=Istrm1,Iendp2
!>            tl_dZde_p(i,j,k2)=0.5_r8*(tl_VFe(i-1,j)+                  &
!>   &                                  tl_VFe(i  ,j))
!>
              adfac=0.5_r8*ad_dZde_p(i,j,k2)
              ad_VFe(i-1,j)=ad_VFe(i-1,j)+adfac
              ad_VFe(i  ,j)=ad_VFe(i  ,j)+adfac
              ad_dZde_p(i,j,k2)=0.0_r8
!>            tl_dZdx_p(i,j,k2)=0.5_r8*(tl_UFx(i,j-1)+                  &
!>   &                                  tl_UFx(i,j  ))
!>
              adfac=0.5_r8*ad_dZdx_p(i,j,k2)
              ad_UFx(i,j-1)=ad_UFx(i,j-1)+adfac
              ad_UFx(i,j  )=ad_UFx(i,j  )+adfac
              ad_dZdx_p(i,j,k2)=0.0_r8
            END DO
          END DO
!
          DO j=JstrVm2,Jendp2
            DO i=Istrm2,Iendp2
              cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
!>            tl_VFe(i,j)=cff*(tl_z_r(i,j  ,k+1)-                       &
!>   &                         tl_z_r(i,j-1,k+1))
!>
              adfac=cff*ad_VFe(i,j)
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_VFe(i,j)=0.0_r8
            END DO
          END DO

          DO j=Jstrm2,Jendp2
            DO i=IstrUm2,Iendp2
              cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
!>            tl_UFx(i,j)=cff*(tl_z_r(i  ,j,k+1)-                       &
!>   &                         tl_z_r(i-1,j,k+1))
!>
              adfac=cff*ad_UFx(i,j)
              ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-adfac
              ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
              ad_UFx(i,j)=0.0_r8
            END DO
          END DO
        END IF
!
!  Compute new recursive storage indices.
!
        kt=k2
        k2=k1
        k1=kt
      END DO K_LOOP3
      RETURN
      END SUBROUTINE ad_uv3dmix4_tile
