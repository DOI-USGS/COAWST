      SUBROUTINE ad_t3dmix2 (ng, tile)
!
!svn $Id: ad_t3dmix2_geo.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes adjoint horizontal harmonic mixing of      !
!  tracers along geopotential surfaces.                                !
!                                                                      !
!  BASIC STATE variables needed: diff2, Hz, t, z_r                     !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef TS_MIX_CLIMA
      USE mod_clima
#endif
#ifdef DIAGNOSTICS_TS
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
      CALL wclock_on (ng, iADM, 25)
#endif
      CALL ad_t3dmix2_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nstp(ng), nnew(ng),               &
#ifdef MASKING
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
#ifdef WET_DRY_NOT_YET
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
#endif
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % ad_z_r,                          &
#ifdef DIFF_3DCOEF
     &                      MIXING(ng) % diff3d_r,                      &
#else
     &                      MIXING(ng) % diff2,                         &
#endif
#ifdef TS_MIX_CLIMA
     &                      CLIMA(ng) % tclm,                           &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                      DIAGS(ng) % DiaTwrk,                        &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t)
#ifdef PROFILE
      CALL wclock_off (ng, iADM, 25)
#endif
      RETURN
      END SUBROUTINE ad_t3dmix2
!
!***********************************************************************
      SUBROUTINE ad_t3dmix2_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
#ifdef WET_DRY_NOT_YET
     &                            umask_wet, vmask_wet,                 &
#endif
     &                            om_v, on_u, pm, pn,                   &
     &                            Hz, ad_Hz,                            &
     &                            z_r, ad_z_r,                          &
#ifdef DIFF_3DCOEF
     &                            diff3d_r,                             &
#else
     &                            diff2,                                &
#endif
#ifdef TS_MIX_CLIMA
     &                            tclm,                                 &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                            DiaTwrk,                              &
#endif
     &                            t, ad_t)
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
      integer, intent(in) :: nrhs, nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
# ifdef DIFF_3DCOEF
      real(r8), intent(in) :: diff3d_r(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
# endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef DIFF_3DCOEF
      real(r8), intent(in) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
# endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k, kk, kt, k1, k1b, k2, k2b

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3, ad_cff4
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_FS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dTdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dZde

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_cff4=0.0_r8

      ad_FE(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS,JminS:JmaxS)=0.0_r8

      ad_FS(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8

      ad_dTdz(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dTdx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dTde(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dZdx(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
      ad_dZde(IminS:ImaxS,JminS:JmaxS,1:2)=0.0_r8
!
!----------------------------------------------------------------------
!  Compute horizontal harmonic diffusion along geopotential surfaces.
!----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!        dTdx,dTde(:,:,k1) k     rho-points
!        dTdx,dTde(:,:,k2) k+1   rho-points
!          FS,dTdz(:,:,k1) k-1/2   W-points
!          FS,dTdz(:,:,k2) k+1/2   W-points

#ifdef TS_MIX_STABILITY
!
!  In order to increase stability, the biharmonic operator is applied
!  as: 3/4 t(:,:,:,nrhs,:) + 1/4 t(:,:,:,nstp,:).
#endif
!
!  Compute adjoint of starting values of k1 and k2.
!
      T_LOOP : DO itrc=1,NT(ng)
        k1=2
        k2=1
        DO k=0,N(ng)
!!
!!  Note: The following code is equivalent to
!!
!!        kt=k1
!!        k1=k2
!!        k2=kt
!!
!!  We use the adjoint of above code.
!!
          k1=k2
          k2=3-k1
        END DO
        K_LOOP : DO k=N(ng),0,-1
!
!  Compute required BASIC STATE fields. Need to look forward in
!  recursive kk index.
!
          k2b=1
          DO kk=0,k
            k1b=k2b
            k2b=3-k1b
!
!  Compute components of the rotated tracer flux (T m3/s) along
!  geopotential surfaces (required BASIC STATE fields).
!
            IF (kk.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend+1
                  cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                  cff=cff*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                  cff=cff*umask_wet(i,j)
#endif
                  dZdx(i,j,k2b)=cff*(z_r(i  ,j,kk+1)-                   &
     &                               z_r(i-1,j,kk+1))
#if defined TS_MIX_STABILITY
                  dTdx(i,j,k2b)=cff*(0.75_r8*(t(i  ,j,kk+1,nrhs,itrc)-  &
     &                                        t(i-1,j,kk+1,nrhs,itrc))+ &
     &                               0.25_r8*(t(i  ,j,kk+1,nstp,itrc)-  &
     &                                        t(i-1,j,kk+1,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                  IF (LtracerCLM(itrc,ng)) THEN
                    dTdx(i,j,k2b)=cff*((t(i  ,j,kk+1,nrhs,itrc)-        &
     &                                  tclm(i  ,j,kk+1,itrc))-         &
     &                                 (t(i-1,j,kk+1,nrhs,itrc)-        &
     &                                  tclm(i-1,j,kk+1,itrc)))
                  ELSE
                    dTdx(i,j,k2b)=cff*(t(i  ,j,kk+1,nrhs,itrc)-         &
     &                                 t(i-1,j,kk+1,nrhs,itrc))
                  END IF
#else
                  dTdx(i,j,k2b)=cff*(t(i  ,j,kk+1,nrhs,itrc)-           &
     &                               t(i-1,j,kk+1,nrhs,itrc))
#endif
                END DO
              END DO
              IF (kk.eq.0) THEN
                DO j=Jstr,Jend
                  DO i=Istr,Iend+1
                    dZdx(i,j,k1b)=0.0_r8
                    dTdx(i,j,k1b)=0.0_r8
                  END DO
                END DO
              END IF
              DO j=Jstr,Jend+1
                DO i=Istr,Iend
                  cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                  cff=cff*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                  cff=cff*vmask_wet(i,j)
#endif
                  dZde(i,j,k2b)=cff*(z_r(i,j  ,kk+1)-                   &
     &                               z_r(i,j-1,kk+1))
#if defined TS_MIX_STABILITY
                  dTde(i,j,k2b)=cff*(0.75_r8*(t(i,j  ,kk+1,nrhs,itrc)-  &
     &                                        t(i,j-1,kk+1,nrhs,itrc))+ &
     &                               0.25_r8*(t(i,j  ,kk+1,nstp,itrc)-  &
     &                                        t(i,j-1,kk+1,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                  IF (LtracerCLM(itrc,ng)) THEN
                    dTde(i,j,k2b)=cff*((t(i,j  ,kk+1,nrhs,itrc)-        &
     &                                  tclm(i,j  ,kk+1,itrc))-         &
     &                                 (t(i,j-1,kk+1,nrhs,itrc)-        &
     &                                  tclm(i,j-1,kk+1,itrc)))
                  ELSE
                    dTde(i,j,k2b)=cff*(t(i,j  ,kk+1,nrhs,itrc)-         &
     &                                 t(i,j-1,kk+1,nrhs,itrc))
                  END IF
#else
                  dTde(i,j,k2b)=cff*(t(i,j  ,kk+1,nrhs,itrc)-           &
     &                               t(i,j-1,kk+1,nrhs,itrc))
#endif
                END DO
              END DO
              IF (kk.eq.0) THEN
                DO j=Jstr,Jend+1
                  DO i=Istr,Iend
                    dZde(i,j,k1b)=0.0_r8
                    dTde(i,j,k1b)=0.0_r8
                  END DO
                END DO
              END IF
            END IF
            IF ((kk.eq.0).or.(kk.eq.N(ng))) THEN
              DO j=Jstr-1,Jend+1
                DO i=Istr-1,Iend+1
                  dTdz(i,j,k2b)=0.0_r8
                END DO
              END DO
              IF (kk.eq.0) THEN
                DO j=Jstr-1,Jend+1
                  DO i=Istr-1,Iend+1
                    dTdz(i,j,k1b)=0.0_r8
                  END DO
                END DO
              END IF
            ELSE
              DO j=Jstr-1,Jend+1
                DO i=Istr-1,Iend+1
                  cff=1.0_r8/(z_r(i,j,kk+1)-z_r(i,j,kk))
#if defined TS_MIX_STABILITY
                  dTdz(i,j,k2b)=cff*(0.75_r8*(t(i,j,kk+1,nrhs,itrc)-    &
     &                                        t(i,j,kk  ,nrhs,itrc))+   &
     &                               0.25_r8*(t(i,j,kk+1,nstp,itrc)-    &
     &                                        t(i,j,kk  ,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                  IF (LtracerCLM(itrc,ng)) THEN
                    dTdz(i,j,k2b)=cff*((t(i,j,kk+1,nrhs,itrc)-          &
     &                                  tclm(i,j,kk+1,itrc))-           &
     &                                 (t(i,j,kk  ,nrhs,itrc)-          &
     &                                  tclm(i,j,kk  ,itrc)))
                  ELSE
                    dTdz(i,j,k2b)=cff*(t(i,j,kk+1,nrhs,itrc)-           &
     &                                 t(i,j,kk  ,nrhs,itrc))
                  END IF
#else
                  dTdz(i,j,k2b)=cff*(t(i,j,kk+1,nrhs,itrc)-             &
     &                               t(i,j,kk  ,nrhs,itrc))
#endif
                END DO
              END DO
            END IF
          END DO
!
          IF (k.gt.0) THEN
!
!  Time-step adjoint harmonic, geopotential diffusion term.
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
#ifdef DIAGNOSTICS_TS
!!              DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
!>              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
!>
                ad_cff=ad_cff+ad_t(i,j,k,nnew,itrc)
!>              tl_cff=dt(ng)*pm(i,j)*pn(i,j)*                          &
!>   &                        (tl_FX(i+1,j)-tl_FX(i,j)+                 &
!>   &                         tl_FE(i,j+1)-tl_FE(i,j))+                &
!>   &                 dt(ng)*(tl_FS(i,j,k2)-tl_FS(i,j,k1))
!>
                adfac=dt(ng)*ad_cff
                adfac1=adfac*pm(i,j)*pn(i,j)
                ad_FS(i,j,k2)=ad_FS(i,j,k2)+adfac
                ad_FS(i,j,k1)=ad_FS(i,j,k1)-adfac
                ad_FE(i,j  )=ad_FE(i,j  )-adfac1
                ad_FE(i,j+1)=ad_FE(i,j+1)+adfac1
                ad_FX(i  ,j)=ad_FX(i  ,j)-adfac1
                ad_FX(i+1,j)=ad_FX(i+1,j)+adfac1
                ad_cff=0.0_r8
              END DO
            END DO
!
!  Compute components of the rotated tracer flux (T m4/s) along
!  geopotential surfaces.
!
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
#ifdef DIFF_3DCOEF
                  cff=0.5_r8*diff3d_r(i,j,k)
#else
                  cff=0.5_r8*diff2(i,j,itrc)
#endif
                  cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde(i,j+1,k1),0.0_r8)
!>                tl_FS(i,j,k2)=tl_FS(i,j,k2)+                          &
!>   &                          cff*                                    &
!>   &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
!>   &                                    dTde(i,j  ,k1))+              &
!>   &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
!>   &                                    dTde(i,j+1,k2))+              &
!>   &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
!>   &                                    dTde(i,j  ,k2))+              &
!>   &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
!>   &                                    dTde(i,j+1,k1))+              &
!>   &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
!>   &                                 cff1*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTde(i,j  ,k1))+              &
!>   &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
!>   &                                 cff2*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTde(i,j+1,k2))+              &
!>   &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
!>   &                                 cff3*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTde(i,j  ,k2))+              &
!>   &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
!>   &                                 cff4*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTde(i,j+1,k1)))
!>
                  adfac=cff*ad_FS(i,j,k2)
                  ad_cff1=ad_cff1+                                      &
     &                    (2.0_r8*cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))*    &
     &                    adfac
                  ad_cff2=ad_cff2+                                      &
     &                    (2.0_r8*cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))*    &
     &                    adfac
                  ad_cff3=ad_cff3+                                      &
     &                    (2.0_r8*cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))*    &
     &                    adfac
                  ad_cff4=ad_cff4+                                      &
     &                    (2.0_r8*cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))*    &
     &                    adfac
                  ad_dTdz(i,j,k2)=ad_dTdz(i,j,k2)+                      &
     &                            (cff1*cff1+                           &
     &                             cff2*cff2+                           &
     &                             cff3*cff3+                           &
     &                             cff4*cff4)*adfac
                  ad_dTde(i,j  ,k1)=ad_dTde(i,j  ,k1)-cff1*adfac
                  ad_dTde(i,j+1,k2)=ad_dTde(i,j+1,k2)-cff2*adfac
                  ad_dTde(i,j  ,k2)=ad_dTde(i,j  ,k2)-cff3*adfac
                  ad_dTde(i,j+1,k1)=ad_dTde(i,j+1,k1)-cff4*adfac
!>                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde(i,j+1,k1)))*        &
!>   &                    tl_dZde(i,j+1,k1)
!>
                  ad_dZde(i,j+1,k1)=ad_dZde(i,j+1,k1)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                            dZde(i,j+1,k1)))*     &
     &                              ad_cff4
                  ad_cff4=0.0_r8
!>                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde(i,j  ,k2)))*        &
!>   &                    tl_dZde(i,j  ,k2)
!>
                  ad_dZde(i,j  ,k2)=ad_dZde(i,j  ,k2)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                            dZde(i,j  ,k2)))*     &
     &                              ad_cff3
                  ad_cff3=0.0_r8
!>                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde(i,j+1,k2)))*        &
!>   &                    tl_dZde(i,j+1,k2)
!>
                  ad_dZde(i,j+1,k2)=ad_dZde(i,j+1,k2)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                           -dZde(i,j+1,k2)))*     &
     &                              ad_cff2
                  ad_cff2=0.0_r8
!>                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde(i,j  ,k1)))*        &
!>   &                    tl_dZde(i,j  ,k1)
!>
                  ad_dZde(i,j  ,k1)=ad_dZde(i,j  ,k1)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                           -dZde(i,j  ,k1)))*     &
     &                               ad_cff1
                  ad_cff1=0.0_r8
!
                  cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
!>                tl_FS(i,j,k2)=cff*                                    &
!>   &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
!>   &                                    dTdx(i  ,j,k1))+              &
!>   &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
!>   &                                    dTdx(i+1,j,k2))+              &
!>   &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
!>   &                                    dTdx(i  ,j,k2))+              &
!>   &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
!>   &                                    dTdx(i+1,j,k1))+              &
!>   &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
!>   &                                 cff1*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTdx(i  ,j,k1))+              &
!>   &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
!>   &                                 cff2*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTdx(i+1,j,k2))+              &
!>   &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
!>   &                                 cff3*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTdx(i  ,j,k2))+              &
!>   &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
!>   &                                 cff4*tl_dTdz(i,j,k2)-            &
!>   &                                 tl_dTdx(i+1,j,k1)))
!>
                  ad_cff1=ad_cff1+                                      &
     &                    (2.0_r8*cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))*    &
     &                    adfac
                  ad_cff2=ad_cff2+                                      &
     &                    (2.0_r8*cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))*    &
     &                    adfac
                  ad_cff3=ad_cff3+                                      &
     &                    (2.0_r8*cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))*    &
     &                    adfac
                  ad_cff4=ad_cff4+                                      &
     &                    (2.0_r8*cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))*    &
     &                    adfac
                  ad_dTdz(i,j,k2)=ad_dTdz(i,j,k2)+                      &
     &                            (cff1*cff1+                           &
     &                             cff2*cff2+                           &
     &                             cff3*cff3+                           &
     &                             cff4*cff4)*adfac
                  ad_dTdx(i  ,j,k1)=ad_dTdx(i  ,j,k1)-cff1*adfac
                  ad_dTdx(i+1,j,k2)=ad_dTdx(i+1,j,k2)-cff2*adfac
                  ad_dTdx(i  ,j,k2)=ad_dTdx(i  ,j,k2)-cff3*adfac
                  ad_dTdx(i+1,j,k1)=ad_dTdx(i+1,j,k1)-cff4*adfac
                  ad_FS(i,j,k2)=0.0_r8
!>                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx(i+1,j,k1)))*        &
!>   &                    tl_dZdx(i+1,j,k1)
!>
                  ad_dZdx(i+1,j,k1)=ad_dZdx(i+1,j,k1)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                            dZdx(i+1,j,k1)))*     &
     &                              ad_cff4
                  ad_cff4=0.0_r8
!>                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx(i  ,j,k2)))*        &
!>   &                    tl_dZdx(i  ,j,k2)
!>
                  ad_dZdx(i  ,j,k2)=ad_dZdx(i  ,j,k2)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                            dZdx(i  ,j,k2)))*     &
     &                              ad_cff3
                  ad_cff3=0.0_r8
!>                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx(i+1,j,k2)))*        &
!>   &                    tl_dZdx(i+1,j,k2)
!>
                  ad_dZdx(i+1,j,k2)=ad_dZdx(i+1,j,k2)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                           -dZdx(i+1,j,k2)))*     &
     &                              ad_cff2
                  ad_cff2=0.0_r8
!>                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx(i  ,j,k1)))*        &
!>   &                    tl_dZdx(i  ,j,k1)
!>
                  ad_dZdx(i  ,j,k1)=ad_dZdx(i  ,j,k1)+                  &
     &                              (0.5_r8+SIGN(0.5_r8,                &
     &                                           -dZdx(i  ,j,k1)))*     &
     &                              ad_cff1
                  ad_cff1=0.0_r8
                END DO
              END DO
            END IF
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
#ifdef DIFF_3DCOEF
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*        &
     &              om_v(i,j)
#else
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*        &
     &              om_v(i,j)
#endif
!>              tl_FE(i,j)=cff*                                         &
!>   &                     (((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*             &
!>   &                       (dTde(i,j,k1)-                             &
!>   &                        0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*         &
!>   &                                   (dTdz(i,j-1,k1)+               &
!>   &                                    dTdz(i,j  ,k2))+              &
!>   &                                MAX(dZde(i,j,k1),0.0_r8)*         &
!>   &                                   (dTdz(i,j-1,k2)+               &
!>   &                                    dTdz(i,j  ,k1)))))+           &
!>   &                      ((Hz(i,j,k)+Hz(i,j-1,k))*                   &
!>   &                       (tl_dTde(i,j,k1)-                          &
!>   &                        0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*         &
!>   &                                   (tl_dTdz(i,j-1,k1)+            &
!>   &                                    tl_dTdz(i,j  ,k2))+           &
!>   &                                MAX(dZde(i,j,k1),0.0_r8)*         &
!>   &                                   (tl_dTdz(i,j-1,k2)+            &
!>   &                                    tl_dTdz(i,j  ,k1)))-          &
!>   &                        0.5_r8*((0.5_r8+                          &
!>   &                                 SIGN(0.5_r8,-dZde(i,j,k1)))*     &
!>   &                                tl_dZde(i,j,k1)*                  &
!>   &                                (dTdz(i,j-1,k1)+dTdz(i,j,k2))+    &
!>   &                                (0.5_r8+                          &
!>   &                                 SIGN(0.5_r8, dZde(i,j,k1)))*     &
!>   &                                tl_dZde(i,j,k1)*                  &
!>   &                                (dTdz(i,j-1,k2)+dTdz(i,j,k1))))))
!>
                adfac=cff*ad_FE(i,j)
                adfac1=adfac*(dTde(i,j,k1)-                             &
     &                        0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*         &
     &                                   (dTdz(i,j-1,k1)+               &
     &                                    dTdz(i,j  ,k2))+              &
     &                                MAX(dZde(i,j,k1),0.0_r8)*         &
     &                                   (dTdz(i,j-1,k2)+               &
     &                                    dTdz(i,j  ,k1))))
                adfac2=adfac*(Hz(i,j,k)+Hz(i,j-1,k))
                adfac3=adfac2*0.5_r8*MIN(dZde(i,j,k1),0.0_r8)
                adfac4=adfac2*0.5_r8*MAX(dZde(i,j,k1),0.0_r8)
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_dTde(i,j,k1)=ad_dTde(i,j,k1)+adfac2
                ad_dTdz(i,j-1,k1)=ad_dTdz(i,j-1,k1)-adfac3
                ad_dTdz(i,j  ,k2)=ad_dTdz(i,j  ,k2)-adfac3
                ad_dTdz(i,j-1,k2)=ad_dTdz(i,j-1,k2)-adfac4
                ad_dTdz(i,j  ,k1)=ad_dTdz(i,j  ,k1)-adfac4
                ad_dZde(i,j,k1)=ad_dZde(i,j,k1)-                        &
     &                          adfac2*0.5_r8*                          &
     &                          ((0.5_r8+SIGN(0.5_r8,-dZde(i,j,k1)))*   &
     &                           (dTdz(i,j-1,k1)+dTdz(i,j,k2))+         &
     &                           (0.5_r8+SIGN(0.5_r8, dZde(i,j,k1)))*   &
     &                           (dTdz(i,j-1,k2)+dTdz(i,j,k1)))
                ad_FE(i,j)=0.0_r8
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
#ifdef DIFF_3DCOEF
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*        &
     &              on_u(i,j)
#else
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*        &
     &              on_u(i,j)
#endif
!>              tl_FX(i,j)=cff*                                         &
!>   &                     (((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*             &
!>   &                       (dTdx(i,j,k1)-                             &
!>   &                        0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*         &
!>   &                                   (dTdz(i-1,j,k1)+               &
!>   &                                    dTdz(i  ,j,k2))+              &
!>   &                                MAX(dZdx(i,j,k1),0.0_r8)*         &
!>   &                                   (dTdz(i-1,j,k2)+               &
!>   &                                    dTdz(i  ,j,k1)))))+           &
!>   &                      ((Hz(i,j,k)+Hz(i-1,j,k))*                   &
!>   &                       (tl_dTdx(i,j,k1)-                          &
!>   &                        0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*         &
!>   &                                   (tl_dTdz(i-1,j,k1)+            &
!>   &                                    tl_dTdz(i  ,j,k2))+           &
!>   &                                MAX(dZdx(i,j,k1),0.0_r8)*         &
!>   &                                   (tl_dTdz(i-1,j,k2)+            &
!>   &                                    tl_dTdz(i  ,j,k1)))-          &
!>   &                        0.5_r8*((0.5_r8+                          &
!>   &                                 SIGN(0.5_r8,-dZdx(i,j,k1)))*     &
!>   &                                tl_dZdx(i,j,k1)*                  &
!>   &                                (dTdz(i-1,j,k1)+dTdz(i,j,k2))+    &
!>   &                                (0.5_r8+                          &
!>   &                                 SIGN(0.5_r8, dZdx(i,j,k1)))*     &
!>   &                                tl_dZdx(i,j,k1)*                  &
!>   &                                (dTdz(i-1,j,k2)+dTdz(i,j,k1))))))
!>
                adfac=cff*ad_FX(i,j)
                adfac1=adfac*(dTdx(i,j,k1)-                             &
     &                        0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*         &
     &                                   (dTdz(i-1,j,k1)+               &
     &                                    dTdz(i  ,j,k2))+              &
     &                                MAX(dZdx(i,j,k1),0.0_r8)*         &
     &                                   (dTdz(i-1,j,k2)+               &
     &                                    dTdz(i  ,j,k1))))
                adfac2=adfac*(Hz(i,j,k)+Hz(i-1,j,k))
                adfac3=adfac2*0.5_r8*MIN(dZdx(i,j,k1),0.0_r8)
                adfac4=adfac2*0.5_r8*MAX(dZdx(i,j,k1),0.0_r8)
                ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
                ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
                ad_dTdx(i,j,k1)=ad_dTdx(i,j,k1)+adfac2
                ad_dTdz(i-1,j,k1)=ad_dTdz(i-1,j,k1)-adfac3
                ad_dTdz(i  ,j,k2)=ad_dTdz(i  ,j,k2)-adfac3
                ad_dTdz(i-1,j,k2)=ad_dTdz(i-1,j,k2)-adfac4
                ad_dTdz(i  ,j,k1)=ad_dTdz(i  ,j,k1)-adfac4
                ad_dZdx(i,j,k1)=ad_dZdx(i,j,k1)-                        &
     &                          adfac2*0.5_r8*                          &
     &                          ((0.5_r8+SIGN(0.5_r8,-dZdx(i,j,k1)))*   &
     &                           (dTdz(i-1,j,k1)+dTdz(i,j,k2))+         &
     &                           (0.5_r8+SIGN(0.5_r8, dZdx(i,j,k1)))*   &
     &                           (dTdz(i-1,j,k2)+dTdz(i,j,k1)))
                ad_FX(i,j)=0.0_r8
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
!>              tl_FS(i,j,k2)=0.0_r8
!>
                ad_FS(i,j,k2)=0.0_r8
!>              tl_dTdz(i,j,k2)=0.0_r8
!>
                ad_dTdz(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
#if defined TS_MIX_STABILITY
!>              tl_dTdz(i,j,k2)=tl_cff*(0.75_r8*(t(i,j,k+1,nrhs,itrc)-  &
!>   &                                           t(i,j,k  ,nrhs,itrc))+ &
!>   &                                  0.25_r8*(t(i,j,k+1,nstp,itrc)-  &
!>   &                                           t(i,j,k  ,nstp,itrc)))+&
!>   &                          cff*(0.75_r8*(tl_t(i,j,k+1,nrhs,itrc)-  &
!>   &                                        tl_t(i,j,k  ,nrhs,itrc))+ &
!>   &                               0.25_r8*(tl_t(i,j,k+1,nstp,itrc)-  &
!>   &                                        tl_t(i,j,k  ,nstp,itrc)))
!>
                adfac=cff*ad_dTdz(i,j,k2)
                adfac1=adfac*0.75_r8
                adfac2=adfac*0.25_r8
                ad_t(i,j,k  ,nrhs,itrc)=ad_t(i,j,k  ,nrhs,itrc)-adfac1
                ad_t(i,j,k+1,nrhs,itrc)=ad_t(i,j,k+1,nrhs,itrc)+adfac1
                ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)-adfac2
                ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+adfac2
                ad_cff=ad_cff+(0.75_r8*(t(i,j,k+1,nrhs,itrc)-           &
     &                                  t(i,j,k  ,nrhs,itrc))+          &
     &                         0.25_r8*(t(i,j,k+1,nstp,itrc)-           &
     &                                  t(i,j,k  ,nstp,itrc)))*         &
     &                        ad_dTdz(i,j,k2)
                ad_dTdz(i,j,k2)=0.0_r8
#elif defined TS_MIX_CLIMA
                IF (LtracerCLM(itrc,ng)) THEN
!>                tl_dTdz(i,j,k2)=tl_cff*((t(i,j,k+1,nrhs,itrc)-        &
!>   &                                     tclm(i,j,k+1,itrc))-         &
!>   &                                    (t(i,j,k  ,nrhs,itrc)-        &
!>   &                                     tclm(i,j,k  ,itrc)))+        &
!>   &                            cff*(tl_t(i,j,k+1,nrhs,itrc)-         &
!>   &                                 tl_t(i,j,k  ,nrhs,itrc))
!>
                  adfac=cff*ad_dTdz(i,j,k2)
                  ad_t(i,j,k  ,nrhs,itrc)=ad_t(i,j,k  ,nrhs,itrc)-adfac
                  ad_t(i,j,k+1,nrhs,itrc)=ad_t(i,j,k+1,nrhs,itrc)+adfac
                  ad_cff=ad_cff+((t(i,j,k+1,nrhs,itrc)-                 &
     &                            tclm(i,j,k+1,itrc))-                  &
     &                           (t(i,j,k  ,nrhs,itrc)-                 &
     &                            tclm(i,j,k  ,itrc)))*ad_dTdz(i,j,k2)
                  ad_dTdz(i,j,k2)=0.0_r8
                ELSE
!>                tl_dTdz(i,j,k2)=tl_cff*(t(i,j,k+1,nrhs,itrc)-         &
!>   &                                    t(i,j,k  ,nrhs,itrc))+        &
!>   &                            cff*(tl_t(i,j,k+1,nrhs,itrc)-         &
!>   &                                 tl_t(i,j,k  ,nrhs,itrc))
!>
                  adfac=cff*ad_dTdz(i,j,k2)
                  ad_t(i,j,k  ,nrhs,itrc)=ad_t(i,j,k  ,nrhs,itrc)-adfac
                  ad_t(i,j,k+1,nrhs,itrc)=ad_t(i,j,k+1,nrhs,itrc)+adfac
                  ad_cff=ad_cff+(t(i,j,k+1,nrhs,itrc)-                  &
     &                           t(i,j,k  ,nrhs,itrc))*ad_dTdz(i,j,k2)
                  ad_dTdz(i,j,k2)=0.0_r8
                END IF
#else
!>              tl_dTdz(i,j,k2)=tl_cff*(t(i,j,k+1,nrhs,itrc)-           &
!>   &                                  t(i,j,k  ,nrhs,itrc))+          &
!>   &                          cff*(tl_t(i,j,k+1,nrhs,itrc)-           &
!>   &                               tl_t(i,j,k  ,nrhs,itrc))
!>
                adfac=cff*ad_dTdz(i,j,k2)
                ad_t(i,j,k  ,nrhs,itrc)=ad_t(i,j,k  ,nrhs,itrc)-adfac
                ad_t(i,j,k+1,nrhs,itrc)=ad_t(i,j,k+1,nrhs,itrc)+adfac
                ad_cff=ad_cff+(t(i,j,k+1,nrhs,itrc)-                    &
     &                         t(i,j,k  ,nrhs,itrc))*ad_dTdz(i,j,k2)
                ad_dTdz(i,j,k2)=0.0_r8
#endif
!>              tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
!>
                adfac=cff*cff*ad_cff
                ad_z_r(i,j,k  )=ad_z_r(i,j,k  )+adfac
                ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)-adfac
                ad_cff=0.0_r8
              END DO
            END DO
          END IF
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*vmask_wet(i,j)
#endif
#if defined TS_MIX_STABILITY
!>              tl_dTde(i,j,k2)=cff*                                    &
!>   &                          (0.75_r8*(tl_t(i,j  ,k+1,nrhs,itrc)-    &
!>   &                                    tl_t(i,j-1,k+1,nrhs,itrc))+   &
!>   &                           0.25_r8*(tl_t(i,j  ,k+1,nstp,itrc)-    &
!>   &                                    tl_t(i,j-1,k+1,nstp,itrc)))
!>
                adfac=cff*ad_dTde(i,j,k2)
                adfac1=adfac*0.75_r8
                adfac2=adfac*0.25_r8
                ad_t(i,j-1,k+1,nrhs,itrc)=ad_t(i,j-1,k+1,nrhs,itrc)-    &
     &                                    adfac1
                ad_t(i,j  ,k+1,nrhs,itrc)=ad_t(i,j  ,k+1,nrhs,itrc)+    &
     &                                    adfac1
                ad_t(i,j-1,k+1,nstp,itrc)=ad_t(i,j-1,k+1,nstp,itrc)-    &
     &                                    adfac2
                ad_t(i,j  ,k+1,nstp,itrc)=ad_t(i,j  ,k+1,nstp,itrc)+    &
     &                                    adfac2
                ad_dTde(i,j,k2)=0.0_r8
#elif defined TS_MIX_CLIMA
!>              tl_dTde(i,j,k2)=cff*(tl_t(i,j  ,k+1,nrhs,itrc)-         &
!>   &                               tl_t(i,j-1,k+1,nrhs,itrc))
!>
                adfac=cff*ad_dTde(i,j,k2)
                ad_t(i,j-1,k+1,nrhs,itrc)=ad_t(i,j-1,k+1,nrhs,itrc)-    &
     &                                    adfac
                ad_t(i,j  ,k+1,nrhs,itrc)=ad_t(i,j  ,k+1,nrhs,itrc)+    &
     &                                    adfac
                ad_dTde(i,j,k2)=0.0_r8
#else
!>              tl_dTde(i,j,k2)=cff*(tl_t(i,j  ,k+1,nrhs,itrc)-         &
!>   &                               tl_t(i,j-1,k+1,nrhs,itrc))
!>
                adfac=cff*ad_dTde(i,j,k2)
                ad_t(i,j-1,k+1,nrhs,itrc)=ad_t(i,j-1,k+1,nrhs,itrc)-    &
     &                                    adfac
                ad_t(i,j  ,k+1,nrhs,itrc)=ad_t(i,j  ,k+1,nrhs,itrc)+    &
     &                                    adfac
                ad_dTde(i,j,k2)=0.0_r8
#endif
!>              tl_dZde(i,j,k2)=cff*(tl_z_r(i,j  ,k+1)-                 &
!>   &                               tl_z_r(i,j-1,k+1))
!>
                adfac=cff*ad_dZde(i,j,k2)
                ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-adfac
                ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
                ad_dZde(i,j,k2)=0.0_r8
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*umask_wet(i,j)
#endif
#if defined TS_MIX_STABILITY
!>              tl_dTdx(i,j,k2)=cff*                                    &
!>   &                          (0.75_r8*(tl_t(i  ,j,k+1,nrhs,itrc)-    &
!>   &                                    tl_t(i-1,j,k+1,nrhs,itrc))+   &
!>   &                           0.25_r8*(tl_t(i  ,j,k+1,nstp,itrc)-    &
!>   &                                    tl_t(i-1,j,k+1,nstp,itrc)))
!>
                adfac=cff*ad_dTdx(i,j,k2)
                adfac1=adfac*0.75_r8
                adfac2=adfac*0.25_r8
                ad_t(i-1,j,k+1,nrhs,itrc)=ad_t(i-1,j,k+1,nrhs,itrc)-    &
     &                                    adfac1
                ad_t(i  ,j,k+1,nrhs,itrc)=ad_t(i  ,j,k+1,nrhs,itrc)+    &
     &                                    adfac1
                ad_t(i-1,j,k+1,nstp,itrc)=ad_t(i-1,j,k+1,nstp,itrc)-    &
     &                                    adfac1
                ad_t(i  ,j,k+1,nstp,itrc)=ad_t(i  ,j,k+1,nstp,itrc)+    &
     &                                    adfac1
                ad_dTdx(i,j,k2)=0.0_r8
#elif defined TS_MIX_CLIMA
!>              tl_dTdx(i,j,k2)=cff*(tl_t(i  ,j,k+1,nrhs,itrc)-         &
!>   &                               tl_t(i-1,j,k+1,nrhs,itrc))
!>
                adfac=cff*ad_dTdx(i,j,k2)
                ad_t(i-1,j,k+1,nrhs,itrc)=ad_t(i-1,j,k+1,nrhs,itrc)-    &
     &                                    adfac
                ad_t(i  ,j,k+1,nrhs,itrc)=ad_t(i  ,j,k+1,nrhs,itrc)+    &
     &                                    adfac
                ad_dTdx(i,j,k2)=0.0_r8
#else
!>              tl_dTdx(i,j,k2)=cff*(tl_t(i  ,j,k+1,nrhs,itrc)-         &
!>   &                               tl_t(i-1,j,k+1,nrhs,itrc))
!>
                adfac=cff*ad_dTdx(i,j,k2)
                ad_t(i-1,j,k+1,nrhs,itrc)=ad_t(i-1,j,k+1,nrhs,itrc)-    &
     &                                    adfac
                ad_t(i  ,j,k+1,nrhs,itrc)=ad_t(i  ,j,k+1,nrhs,itrc)+    &
     &                                    adfac
                ad_dTdx(i,j,k2)=0.0_r8
#endif
!>              tl_dZdx(i,j,k2)=cff*(tl_z_r(i  ,j,k+1)-                 &
!>   &                               tl_z_r(i-1,j,k+1))
!>
                adfac=cff*ad_dZdx(i,j,k2)
                ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-adfac
                ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
                ad_dZdx(i,j,k2)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute new storage recursive indices.
!
          kt=k2
          k2=k1
          k1=kt
        END DO K_LOOP
      END DO T_LOOP

      RETURN
      END SUBROUTINE ad_t3dmix2_tile
