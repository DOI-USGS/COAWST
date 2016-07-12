      SUBROUTINE tl_t3dmix4 (ng, tile)
!
!svn $Id: tl_t3dmix4_geo.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes tangent linear horizontal biharmonic       !
!  mixing of tracers along geopotential surfaces.                      !
!                                                                      !
!  BASIC STATE variables needed: diff4, Hz, t, z_r                     !
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
      CALL wclock_on (ng, iTLM, 28)
#endif
      CALL tl_t3dmix4_tile (ng, tile,                                   &
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
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % tl_z_r,                          &
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
     &                      MIXING(ng) % diff3d_u,                      &
     &                      MIXING(ng) % diff3d_v,                      &
# else
     &                      MIXING(ng) % diff3d_r,                      &
# endif
#else
     &                      MIXING(ng) % diff4,                         &
#endif
#ifdef TS_MIX_CLIMA
     &                      CLIMA(ng) % tclm,                           &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                      DIAGS(ng) % DiaTwrk,                        &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)
#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 28)
#endif

      RETURN
      END SUBROUTINE tl_t3dmix4
!
!***********************************************************************
      SUBROUTINE tl_t3dmix4_tile (ng, tile,                             &
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
     &                            Hz, tl_Hz,                            &
     &                            z_r, tl_z_r,                          &
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
     &                            diff3d_u, diff3d_v,                   &
# else
     &                            diff3d_r,                             &
# endif
#else
     &                            diff4,                                &
#endif
#ifdef TS_MIX_CLIMA
     &                            tclm,                                 &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                            DiaTwrk,                              &
#endif
     &                            t, tl_t)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
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
#  ifdef TS_U3ADV_SPLIT
      real(r8), intent(in) :: diff3d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: diff3d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: diff3d_r(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:,LBj:,:)
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
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
      real(r8), intent(in) :: diff3d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: diff3d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
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
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k, k1, k2

      real(r8) :: cff, cff1, cff2, cff3, cff4, dife, difx
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3, tl_cff4

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: LapT

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_LapT

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: FS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_FS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZdx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute tangent linear horizontal biharmonic diffusion along
!  geopotential surfaces.  The biharmonic operator is computed by
!  applying the harmonic operator twice.
!-----------------------------------------------------------------------
!
!  Set local I- and J-ranges.
!
      IF (EWperiodic(ng)) THEN
        Imin=Istr-1
        Imax=Iend+1
      ELSE
        Imin=MAX(Istr-1,1)
        Imax=MIN(Iend+1,Lm(ng))
      END IF
      IF (NSperiodic(ng)) THEN
        Jmin=Jstr-1
        Jmax=Jend+1
      ELSE
        Jmin=MAX(Jstr-1,1)
        Jmax=MIN(Jend+1,Mm(ng))
      END IF
!
!  Compute horizontal and vertical gradients associated with the
!  first rotated harmonic operator.  Notice the recursive blocking
!  sequence. The vertical placement of the gradients is:
!
!        dTdx,dTde(:,:,k1) k     rho-points
!        dTdx,dTde(:,:,k2) k+1   rho-points
!          FS,dTdz(:,:,k1) k-1/2   W-points
!          FS,dTdz(:,:,k2) k+1/2   W-points
!
#ifdef TS_MIX_STABILITY
!  In order to increase stability, the rotated biharmonic is applied
!  as: 3/4 t(:,:,:,nrhs,:) + 1/4 t(:,:,:,nstp,:).
!
#endif

      T_LOOP : DO itrc=1,NT(ng)
        k2=1
        K_LOOP1 : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax+1
                cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*umask_wet(i,j)
#endif
                dZdx(i,j,k2)=cff*(z_r(i  ,j,k+1)-                       &
     &                            z_r(i-1,j,k+1))
                tl_dZdx(i,j,k2)=cff*(tl_z_r(i  ,j,k+1)-                 &
     &                               tl_z_r(i-1,j,k+1))
#if defined TS_MIX_STABILITY
                dTdx(i,j,k2)=cff*(0.75_r8*(t(i  ,j,k+1,nrhs,itrc)-      &
     &                                     t(i-1,j,k+1,nrhs,itrc))+     &
     &                            0.25_r8*(t(i  ,j,k+1,nstp,itrc)-      &
     &                                     t(i-1,j,k+1,nstp,itrc)))
                tl_dTdx(i,j,k2)=cff*
     &                          (0.75_r8*(tl_t(i  ,j,k+1,nrhs,itrc)-    &
     &                                    tl_t(i-1,j,k+1,nrhs,itrc))+   &
     &                           0.25_r8*(tl_t(i  ,j,k+1,nstp,itrc)-    &
     &                                    tl_t(i-1,j,k+1,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                IF (LtracerCLM(itrc,ng)) THEN
                  dTdx(i,j,k2)=cff*((t(i  ,j,k+1,nrhs,itrc)-            &
     &                               tclm(i  ,j,k+1,itrc))-             &
     &                              (t(i-1,j,k+1,nrhs,itrc)-            &
     &                               tclm(i-1,j,k+1,itrc)))
                ELSE
                  dTdx(i,j,k2)=cff*(t(i  ,j,k+1,nrhs,itrc)-             &
     &                              t(i-1,j,k+1,nrhs,itrc))
                END IF
                tl_dTdx(i,j,k2)=cff*(tl_t(i  ,j,k+1,nrhs,itrc)-         &
     &                               tl_t(i-1,j,k+1,nrhs,itrc))
#else
                dTdx(i,j,k2)=cff*(t(i  ,j,k+1,nrhs,itrc)-               &
     &                            t(i-1,j,k+1,nrhs,itrc))
                tl_dTdx(i,j,k2)=cff*(tl_t(i  ,j,k+1,nrhs,itrc)-         &
     &                               tl_t(i-1,j,k+1,nrhs,itrc))
#endif
              END DO
            END DO
            DO j=Jmin,Jmax+1
              DO i=Imin,Imax
                cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*vmask_wet(i,j)
#endif
                dZde(i,j,k2)=cff*(z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j-1,k+1))
                tl_dZde(i,j,k2)=cff*(tl_z_r(i,j  ,k+1)-                 &
     &                               tl_z_r(i,j-1,k+1))
#if defined TS_MIX_STABILITY
                dTde(i,j,k2)=cff*(0.75_r8*(t(i,j  ,k+1,nrhs,itrc)-      &
     &                                     t(i,j-1,k+1,nrhs,itrc))+     &
     &                            0.25_r8*(t(i,j  ,k+1,nstp,itrc)-      &
     &                                     t(i,j-1,k+1,nstp,itrc)))
                tl_dTde(i,j,k2)=cff*
     &                          (0.75_r8*(tl_t(i,j  ,k+1,nrhs,itrc)-    &
     &                                    tl_t(i,j-1,k+1,nrhs,itrc))+   &
     &                           0.25_r8*(tl_t(i,j  ,k+1,nstp,itrc)-    &
     &                                    tl_t(i,j-1,k+1,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                IF (LtracerCLM(itrc,ng)) THEN
                  dTde(i,j,k2)=cff*((t(i,j  ,k+1,nrhs,itrc)-            &
     &                               tclm(i,j  ,k+1,itrc))-             &
     &                              (t(i,j-1,k+1,nrhs,itrc)-            &
     &                               tclm(i,j-1,k+1,itrc)))
                ELSE
                  dTde(i,j,k2)=cff*(t(i,j  ,k+1,nrhs,itrc)-             &
     &                              t(i,j-1,k+1,nrhs,itrc))
                END IF
                tl_dTde(i,j,k2)=cff*(tl_t(i,j  ,k+1,nrhs,itrc)-         &
     &                               tl_t(i,j-1,k+1,nrhs,itrc))
#else
                dTde(i,j,k2)=cff*(t(i,j  ,k+1,nrhs,itrc)-               &
     &                            t(i,j-1,k+1,nrhs,itrc))
                tl_dTde(i,j,k2)=cff*(tl_t(i,j  ,k+1,nrhs,itrc)-         &
     &                               tl_t(i,j-1,k+1,nrhs,itrc))
#endif
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jmin-1,Jmax+1
              DO i=Imin-1,Imax+1
                dTdz(i,j,k2)=0.0_r8
                tl_dTdz(i,j,k2)=0.0_r8
                FS(i,j,k2)=0.0_r8
                tl_FS(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jmin-1,Jmax+1
              DO i=Imin-1,Imax+1
                cff=1.0_r8/(z_r(i,j,k+1)-                               &
     &                      z_r(i,j,k  ))
                tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-                       &
     &                           tl_z_r(i,j,k  ))
#if defined TS_MIX_STABILITY
                dTdz(i,j,k2)=cff*(0.75_r8*(t(i,j,k+1,nrhs,itrc)-        &
     &                                     t(i,j,k  ,nrhs,itrc))+       &
     &                            0.25_r8*(t(i,j,k+1,nstp,itrc)-        &
     &                                     t(i,j,k  ,nstp,itrc)))
                tl_dTdz(i,j,k2)=tl_cff*(0.75_r8*(t(i,j,k+1,nrhs,itrc)-  &
     &                                           t(i,j,k  ,nrhs,itrc))+ &
     &                                  0.25_r8*(t(i,j,k+1,nstp,itrc)-  &
     &                                           t(i,j,k  ,nstp,itrc)))+&
     &                          cff*(0.75_r8*(tl_t(i,j,k+1,nrhs,itrc)-  &
     &                                        tl_t(i,j,k  ,nrhs,itrc))+ &
     &                               0.25_r8*(tl_t(i,j,k+1,nstp,itrc)-  &
     &                                        tl_t(i,j,k  ,nstp,itrc)))
#elif defined TS_MIX_CLIMA
                IF (LtracerCLM(itrc,ng)) THEN
                  dTdz(i,j,k2)=cff*((t(i,j,k+1,nrhs,itrc)-              &
     &                               tclm(i,j,k+1,itrc))-               &
     &                              (t(i,j,k  ,nrhs,itrc)-              &
     &                               tclm(i,j,k  ,itrc)))
                  tl_dTdz(i,j,k2)=tl_cff*((t(i,j,k+1,nrhs,itrc)-        &
     &                                     tclm(i,j,k+1,itrc))-         &
     &                                    (t(i,j,k  ,nrhs,itrc)-        &
     &                                     tclm(i,j,k  ,itrc)))+        &
     &                            cff*(tl_t(i,j,k+1,nrhs,itrc)-         &
     &                                 tl_t(i,j,k  ,nrhs,itrc))

                ELSE
                  dTdz(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-               &
     &                              t(i,j,k  ,nrhs,itrc))
                  tl_dTdz(i,j,k2)=tl_cff*(t(i,j,k+1,nrhs,itrc)-         &
     &                                    t(i,j,k  ,nrhs,itrc))+        &
     &                            cff*(tl_t(i,j,k+1,nrhs,itrc)-         &
     &                                 tl_t(i,j,k  ,nrhs,itrc))
                END IF
#else
                dTdz(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-                 &
     &                            t(i,j,k  ,nrhs,itrc))
                tl_dTdz(i,j,k2)=tl_cff*(t(i,j,k+1,nrhs,itrc)-           &
     &                                  t(i,j,k  ,nrhs,itrc))+          &
     &                          cff*(tl_t(i,j,k+1,nrhs,itrc)-           &
     &                               tl_t(i,j,k  ,nrhs,itrc))
#endif
              END DO
            END DO
          END IF
          IF (k.gt.0) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax+1
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                cff=0.5_r8*diff3d_u(i,j,k)*on_u(i,j)
# else
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*        &
     &              on_u(i,j)
# endif
#else
                cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*        &
     &              on_u(i,j)
#endif
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (dTdx(i,j,k1)-                                  &
     &                   0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i-1,j,k1)+                    &
     &                               dTdz(i  ,j,k2))+                   &
     &                           MAX(dZdx(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i-1,j,k2)+                    &
     &                               dTdz(i  ,j,k1))))
                tl_FX(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      (dTdx(i,j,k1)-                              &
     &                       0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i-1,j,k1)+                &
     &                                   dTdz(i  ,j,k2))+               &
     &                               MAX(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i-1,j,k2)+                &
     &                                   dTdz(i  ,j,k1))))+             &
     &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      (tl_dTdx(i,j,k1)-                           &
     &                       0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i-1,j,k1)+             &
     &                                   tl_dTdz(i  ,j,k2))+            &
     &                               MAX(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i-1,j,k2)+             &
     &                                   tl_dTdz(i  ,j,k1)))-           &
     &                       0.5_r8*((0.5_r8+                           &
     &                                SIGN(0.5_r8,-dZdx(i,j,k1)))*      &
     &                               tl_dZdx(i,j,k1)*                   &
     &                               (dTdz(i-1,j,k1)+dTdz(i,j,k2))+     &
     &                               (0.5_r8+                           &
     &                                SIGN(0.5_r8, dZdx(i,j,k1)))*      &
     &                               tl_dZdx(i,j,k1)*                   &
     &                               (dTdz(i-1,j,k2)+dTdz(i,j,k1)))))
              END DO
            END DO
            DO j=Jmin,Jmax+1
              DO i=Imin,Imax
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                cff=0.5_r8*diff3d_v(i,j,k)*om_v(i,j)
# else
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*        &
     &              om_v(i,j)
# endif
#else
                cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*        &
     &              om_v(i,j)
#endif
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (dTde(i,j,k1)-                                  &
     &                   0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i,j-1,k1)+                    &
     &                               dTdz(i,j  ,k2))+                   &
     &                           MAX(dZde(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i,j-1,k2)+                    &
     &                               dTdz(i,j  ,k1))))
                tl_FE(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      (dTde(i,j,k1)-                              &
     &                       0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i,j-1,k1)+                &
     &                                   dTdz(i,j  ,k2))+               &
     &                               MAX(dZde(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i,j-1,k2)+                &
     &                                   dTdz(i,j  ,k1))))+             &
     &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      (tl_dTde(i,j,k1)-                           &
     &                       0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i,j-1,k1)+             &
     &                                   tl_dTdz(i,j  ,k2))+            &
     &                               MAX(dZde(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i,j-1,k2)+             &
     &                                   tl_dTdz(i,j  ,k1)))-           &
     &                       0.5_r8*((0.5_r8+                           &
     &                                SIGN(0.5_r8,-dZde(i,j,k1)))*      &
     &                               tl_dZde(i,j,k1)*                   &
     &                               (dTdz(i,j-1,k1)+dTdz(i,j,k2))+     &
     &                               (0.5_r8+                           &
     &                                SIGN(0.5_r8, dZde(i,j,k1)))*      &
     &                               tl_dZde(i,j,k1)*                   &
     &                               (dTdz(i,j-1,k2)+dTdz(i,j,k1)))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jmin,Jmax
                DO i=Imin,Imax
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                  difx=0.125_r8*(diff3d_u(i,j,k  )+diff3d_u(i+1,j,k  )+ &
     &                           diff3d_u(i,j,k+1)+diff3d_u(i+1,j,k+1))
                  dife=0.125_r8*(diff3d_v(i,j,k  )+diff3d_v(i,j+1,k  )+ &
     &                           diff3d_v(i,j,k+1)+diff3d_v(i,j+1,k+1))
# else
                  difx=0.5_r8*diff3d_r(i,j,k)
                  dife=difx
# endif
#else
                  difx=0.5_r8*diff4(i,j,itrc)
                  dife=difx
#endif
                  cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx(i  ,j,k1)))*        &
     &                    tl_dZdx(i  ,j,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx(i+1,j,k2)))*        &
     &                    tl_dZdx(i+1,j,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx(i  ,j,k2)))*        &
     &                    tl_dZdx(i  ,j,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx(i+1,j,k1)))*        &
     &                    tl_dZdx(i+1,j,k1)
                  FS(i,j,k2)=difx*                                      &
     &                       (cff1*(cff1*dTdz(i,j,k2)-                  &
     &                              dTdx(i  ,j,k1))+                    &
     &                        cff2*(cff2*dTdz(i,j,k2)-                  &
     &                              dTdx(i+1,j,k2))+                    &
     &                        cff3*(cff3*dTdz(i,j,k2)-                  &
     &                              dTdx(i  ,j,k2))+                    &
     &                        cff4*(cff4*dTdz(i,j,k2)-                  &
     &                              dTdx(i+1,j,k1)))
                  tl_FS(i,j,k2)=difx*                                   &
     &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
     &                                    dTdx(i  ,j,k1))+              &
     &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
     &                                    dTdx(i+1,j,k2))+              &
     &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
     &                                    dTdx(i  ,j,k2))+              &
     &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
     &                                    dTdx(i+1,j,k1))+              &
     &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
     &                                 cff1*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i  ,j,k1))+              &
     &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
     &                                 cff2*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i+1,j,k2))+              &
     &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
     &                                 cff3*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i  ,j,k2))+              &
     &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
     &                                 cff4*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i+1,j,k1)))
!
                  cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde(i,j+1,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde(i,j  ,k1)))*        &
     &                    tl_dZde(i,j  ,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde(i,j+1,k2)))*        &
     &                    tl_dZde(i,j+1,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde(i,j  ,k2)))*        &
     &                    tl_dZde(i,j  ,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde(i,j+1,k1)))*        &
     &                    tl_dZde(i,j+1,k1)
                  FS(i,j,k2)=FS(i,j,k2)+                                &
     &                       dife*                                      &
     &                       (cff1*(cff1*dTdz(i,j,k2)-                  &
     &                              dTde(i,j  ,k1))+                    &
     &                        cff2*(cff2*dTdz(i,j,k2)-                  &
     &                              dTde(i,j+1,k2))+                    &
     &                        cff3*(cff3*dTdz(i,j,k2)-                  &
     &                              dTde(i,j  ,k2))+                    &
     &                        cff4*(cff4*dTdz(i,j,k2)-                  &
     &                              dTde(i,j+1,k1)))
                  tl_FS(i,j,k2)=tl_FS(i,j,k2)+                          &
     &                          dife*                                   &
     &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
     &                                    dTde(i,j  ,k1))+              &
     &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
     &                                    dTde(i,j+1,k2))+              &
     &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
     &                                    dTde(i,j  ,k2))+              &
     &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
     &                                    dTde(i,j+1,k1))+              &
     &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
     &                                 cff1*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j  ,k1))+              &
     &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
     &                                 cff2*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j+1,k2))+              &
     &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
     &                                 cff3*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j  ,k2))+              &
     &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
     &                                 cff4*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j+1,k1)))
                END DO
              END DO
            END IF
!
!  Compute first harmonic operator, without mixing coefficient.
!  Multiply by the metrics of the second harmonic operator.  Save
!  into work array "LapT".
!
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                cff=pm(i,j)*pn(i,j)
                cff1=1.0_r8/Hz(i,j,k)
                tl_cff1=-cff1*cff1*tl_Hz(i,j,k)
                LapT(i,j,k)=cff1*(cff*                                  &
     &                            (FX(i+1,j)-FX(i,j)+                   &
     &                             FE(i,j+1)-FE(i,j))+                  &
     &                            (FS(i,j,k2)-FS(i,j,k1)))
                tl_LapT(i,j,k)=tl_cff1*(cff*                            &
     &                                  (FX(i+1,j)-FX(i,j)+             &
     &                                   FE(i,j+1)-FE(i,j))+            &
     &                                  (FS(i,j,k2)-FS(i,j,k1)))+       &
     &                         cff1*(cff*                               &
     &                               (tl_FX(i+1,j)-tl_FX(i,j)+          &
     &                                tl_FE(i,j+1)-tl_FE(i,j))+         &
     &                               (tl_FS(i,j,k2)-tl_FS(i,j,k1)))
              END DO
            END DO
          END IF
        END DO K_LOOP1
!
!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first harmonic operator.
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (tl_LBC(iwest,isTvar(itrc),ng)%closed) THEN
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  LapT(Istr-1,j,k)=0.0_r8
                  tl_LapT(Istr-1,j,k)=0.0_r8
                END DO
              END DO
            ELSE
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  LapT(Istr-1,j,k)=LapT(Istr,j,k)
                  tl_LapT(Istr-1,j,k)=tl_LapT(Istr,j,k)
                END DO
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (tl_LBC(ieast,isTvar(itrc),ng)%closed) THEN
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  LapT(Iend+1,j,k)=0.0_r8
                  tl_LapT(Iend+1,j,k)=0.0_r8
                END DO
              END DO
            ELSE
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  LapT(Iend+1,j,k)=LapT(Iend,j,k)
                  tl_LapT(Iend+1,j,k)=tl_LapT(Iend,j,k)
                END DO
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (tl_LBC(isouth,isTvar(itrc),ng)%closed) THEN
              DO k=1,N(ng)
                DO i=Imin,Imax
                  LapT(i,Jstr-1,k)=0.0_r8
                  tl_LapT(i,Jstr-1,k)=0.0_r8
                END DO
              END DO
            ELSE
              DO k=1,N(ng)
                DO i=Imin,Imax
                  LapT(i,Jstr-1,k)=LapT(i,Jstr,k)
                  tl_LapT(i,Jstr-1,k)=tl_LapT(i,Jstr,k)
                END DO
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (tl_LBC(inorth,isTvar(itrc),ng)%closed) THEN
              DO k=1,N(ng)
                DO i=Imin,Imax
                  LapT(i,Jend+1,k)=0.0_r8
                  tl_LapT(i,Jend+1,k)=0.0_r8
                END DO
              END DO
            ELSE
              DO k=1,N(ng)
                DO i=Imin,Imax
                  LapT(i,Jend+1,k)=LapT(i,Jend,k)
                  tl_LapT(i,Jend+1,k)=tl_LapT(i,Jend,k)
                END DO
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=1,N(ng)
              LapT(Istr-1,Jstr-1,k)=0.5_r8*                             &
     &                              (LapT(Istr  ,Jstr-1,k)+             &
     &                               LapT(Istr-1,Jstr  ,k))
              tl_LapT(Istr-1,Jstr-1,k)=0.5_r8*                          &
     &                                 (tl_LapT(Istr  ,Jstr-1,k)+       &
     &                                  tl_LapT(Istr-1,Jstr  ,k))
            END DO
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=1,N(ng)
              LapT(Iend+1,Jstr-1,k)=0.5_r8*                             &
     &                              (LapT(Iend  ,Jstr-1,k)+             &
     &                               LapT(Iend+1,Jstr  ,k))
              tl_LapT(Iend+1,Jstr-1,k)=0.5_r8*                          &
     &                                 (tl_LapT(Iend  ,Jstr-1,k)+       &
     &                                  tl_LapT(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=1,N(ng)
              LapT(Istr-1,Jend+1,k)=0.5_r8*                             &
     &                              (LapT(Istr  ,Jend+1,k)+             &
     &                               LapT(Istr-1,Jend  ,k))
              tl_LapT(Istr-1,Jend+1,k)=0.5_r8*                          &
     &                                 (tl_LapT(Istr  ,Jend+1,k)+       &
     &                                  tl_LapT(Istr-1,Jend  ,k))
            END DO
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=1,N(ng)
              LapT(Iend+1,Jend+1,k)=0.5_r8*                             &
     &                              (LapT(Iend  ,Jend+1,k)+             &
     &                               LapT(Iend+1,Jend  ,k))
              tl_LapT(Iend+1,Jend+1,k)=0.5_r8*                          &
     &                                 (tl_LapT(Iend  ,Jend+1,k)+       &
     &                                  tl_LapT(Iend+1,Jend  ,k))
            END DO
          END IF
        END IF
!
!  Compute tangent linear horizontal and vertical gradients associated
!  with the second rotated harmonic operator.
!
        k2=1
        K_LOOP2: DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*umask_wet(i,j)
#endif
                dZdx(i,j,k2)=cff*(z_r(i  ,j,k+1)-                       &
     &                            z_r(i-1,j,k+1))
                tl_dZdx(i,j,k2)=cff*(tl_z_r(i  ,j,k+1)-                 &
     &                               tl_z_r(i-1,j,k+1))
                dTdx(i,j,k2)=cff*(LapT(i  ,j,k+1)-                      &
     &                            LapT(i-1,j,k+1))
                tl_dTdx(i,j,k2)=cff*(tl_LapT(i  ,j,k+1)-                &
     &                               tl_LapT(i-1,j,k+1))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
                cff=cff*vmask_wet(i,j)
#endif
                dZde(i,j,k2)=cff*(z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j-1,k+1))
                tl_dZde(i,j,k2)=cff*(tl_z_r(i,j  ,k+1)-                 &
     &                               tl_z_r(i,j-1,k+1))
                dTde(i,j,k2)=cff*(LapT(i,j  ,k+1)-                      &
     &                            LapT(i,j-1,k+1))
                tl_dTde(i,j,k2)=cff*(tl_LapT(i,j  ,k+1)-                &
     &                               tl_LapT(i,j-1,k+1))
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                dTdz(i,j,k2)=0.0_r8
                tl_dTdz(i,j,k2)=0.0_r8
                FS(i,j,k2)=0.0_r8
                tl_FS(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-                               &
     &                      z_r(i,j,k  ))
                tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-                       &
     &                           tl_z_r(i,j,k  ))
                dTdz(i,j,k2)=cff*(LapT(i,j,k+1)-                        &
     &                            LapT(i,j,k  ))
                tl_dTdz(i,j,k2)=tl_cff*(LapT(i,j,k+1)-                  &
     &                                  LapT(i,j,k  ))+                 &
     &                          cff*(tl_LapT(i,j,k+1)-                  &
     &                               tl_LapT(i,j,k  ))
              END DO
            END DO
          END IF
!
!  Compute components of the rotated tracer flux (T m4/s) along
!  geopotential surfaces.
!
          IF (k.gt.0) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                cff=0.5_r8*diff3d_u(i,j,k)*on_u(i,j)
# else
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*        &
     &              on_u(i,j)
# endif
#else
                cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*        &
     &              on_u(i,j)
#endif
!>              FX(i,j)=cff*
!>   &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
!>   &                  (dTdx(i  ,j,k1)-                                &
!>   &                   0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*              &
!>   &                              (dTdz(i-1,j,k1)+                    &
!>   &                               dTdz(i  ,j,k2))+                   &
!>   &                           MAX(dZdx(i,j,k1),0.0_r8)*              &
!>   &                              (dTdz(i-1,j,k2)+                    &
!>   &                               dTdz(i  ,j,k1))))
!>
                tl_FX(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      (dTdx(i,j,k1)-                              &
     &                       0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i-1,j,k1)+                &
     &                                   dTdz(i  ,j,k2))+               &
     &                               MAX(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i-1,j,k2)+                &
     &                                   dTdz(i  ,j,k1))))+             &
     &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      (tl_dTdx(i,j,k1)-                           &
     &                       0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i-1,j,k1)+             &
     &                                   tl_dTdz(i  ,j,k2))+            &
     &                               MAX(dZdx(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i-1,j,k2)+             &
     &                                   tl_dTdz(i  ,j,k1)))-           &
     &                       0.5_r8*((0.5_r8+                           &
     &                                SIGN(0.5_r8,-dZdx(i,j,k1)))*      &
     &                               tl_dZdx(i,j,k1)*                   &
     &                               (dTdz(i-1,j,k1)+dTdz(i,j,k2))+     &
     &                               (0.5_r8+                           &
     &                                SIGN(0.5_r8, dZdx(i,j,k1)))*      &
     &                               tl_dZdx(i,j,k1)*                   &
     &                               (dTdz(i-1,j,k2)+dTdz(i,j,k1)))))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                cff=0.5_r8*diff3d_v(i,j,k)*om_v(i,j)
# else
                cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*        &
     &              om_v(i,j)
# endif
#else
                cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*        &
     &              om_v(i,j)
#endif
!>              FE(i,j)=cff*                                            &
!>   &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
!>   &                  (dTde(i,j,k1)-                                  &
!>   &                   0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*              &
!>   &                              (dTdz(i,j-1,k1)+                    &
!>   &                               dTdz(i,j  ,k2))+                   &
!>   &                           MAX(dZde(i,j,k1),0.0_r8)*              &
!>   &                              (dTdz(i,j-1,k2)+                    &
!>   &                               dTdz(i,j  ,k1))))
!>
                tl_FE(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      (dTde(i,j,k1)-                              &
     &                       0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i,j-1,k1)+                &
     &                                   dTdz(i,j  ,k2))+               &
     &                               MAX(dZde(i,j,k1),0.0_r8)*          &
     &                                  (dTdz(i,j-1,k2)+                &
     &                                   dTdz(i,j  ,k1))))+             &
     &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      (tl_dTde(i,j,k1)-                           &
     &                       0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i,j-1,k1)+             &
     &                                   tl_dTdz(i,j  ,k2))+            &
     &                               MAX(dZde(i,j,k1),0.0_r8)*          &
     &                                  (tl_dTdz(i,j-1,k2)+             &
     &                                   tl_dTdz(i,j  ,k1)))-           &
     &                       0.5_r8*((0.5_r8+                           &
     &                                SIGN(0.5_r8,-dZde(i,j,k1)))*      &
     &                               tl_dZde(i,j,k1)*                   &
     &                               (dTdz(i,j-1,k1)+dTdz(i,j,k2))+     &
     &                               (0.5_r8+                           &
     &                                SIGN(0.5_r8, dZde(i,j,k1)))*      &
     &                               tl_dZde(i,j,k1)*                   &
     &                               (dTdz(i,j-1,k2)+dTdz(i,j,k1)))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
                  difx=0.125_r8*(diff3d_u(i,j,k  )+diff3d_u(i+1,j,k  )+ &
     &                           diff3d_u(i,j,k+1)+diff3d_u(i+1,j,k+1))
                  dife=0.125_r8*(diff3d_v(i,j,k  )+diff3d_v(i,j+1,k  )+ &
     &                           diff3d_v(i,j,k+1)+diff3d_v(i,j+1,k+1))
# else
                  difx=0.5_r8*diff3d_r(i,j,k)
                  dife=difx
# endif
#else
                  difx=0.5_r8*diff4(i,j,itrc)
                  dife=difx
#endif
                  cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx(i  ,j,k1)))*        &
     &                    tl_dZdx(i  ,j,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx(i+1,j,k2)))*        &
     &                    tl_dZdx(i+1,j,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx(i  ,j,k2)))*        &
     &                    tl_dZdx(i  ,j,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx(i+1,j,k1)))*        &
     &                    tl_dZdx(i+1,j,k1)
!>
!>                FS(i,j,k2)=difx*                                      &
!>   &                       (cff1*(cff1*dTdz(i,j,k2)-                  &
!>   &                              dTdx(i  ,j,k1))+                    &
!>   &                        cff2*(cff2*dTdz(i,j,k2)-                  &
!>   &                              dTdx(i+1,j,k2))+                    &
!>   &                        cff3*(cff3*dTdz(i,j,k2)-                  &
!>   &                              dTdx(i  ,j,k2))+                    &
!>   &                        cff4*(cff4*dTdz(i,j,k2)-                  &
!>   &                              dTdx(i+1,j,k1)))
!>
                  tl_FS(i,j,k2)=difx*                                   &
     &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
     &                                    dTdx(i  ,j,k1))+              &
     &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
     &                                    dTdx(i+1,j,k2))+              &
     &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
     &                                    dTdx(i  ,j,k2))+              &
     &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
     &                                    dTdx(i+1,j,k1))+              &
     &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
     &                                 cff1*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i  ,j,k1))+              &
     &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
     &                                 cff2*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i+1,j,k2))+              &
     &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
     &                                 cff3*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i  ,j,k2))+              &
     &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
     &                                 cff4*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTdx(i+1,j,k1)))
!
                  cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde(i,j+1,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde(i,j  ,k1)))*        &
     &                    tl_dZde(i,j  ,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde(i,j+1,k2)))*        &
     &                    tl_dZde(i,j+1,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde(i,j  ,k2)))*        &
     &                    tl_dZde(i,j  ,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde(i,j+1,k1)))*        &
     &                    tl_dZde(i,j+1,k1)
!>                FS(i,j,k2)=FS(i,j,k2)+                                &
!>   &                       dife*                                      &
!>   &                       (cff1*(cff1*dTdz(i,j,k2)-                  &
!>   &                              dTde(i,j  ,k1))+                    &
!>   &                        cff2*(cff2*dTdz(i,j,k2)-                  &
!>   &                              dTde(i,j+1,k2))+                    &
!>   &                        cff3*(cff3*dTdz(i,j,k2)-                  &
!>   &                              dTde(i,j  ,k2))+                    &
!>   &                        cff4*(cff4*dTdz(i,j,k2)-                  &
!>   &                              dTde(i,j+1,k1)))
!>
                  tl_FS(i,j,k2)=tl_FS(i,j,k2)+                          &
     &                          dife*                                   &
     &                          (tl_cff1*(cff1*dTdz(i,j,k2)-            &
     &                                    dTde(i,j  ,k1))+              &
     &                           tl_cff2*(cff2*dTdz(i,j,k2)-            &
     &                                    dTde(i,j+1,k2))+              &
     &                           tl_cff3*(cff3*dTdz(i,j,k2)-            &
     &                                    dTde(i,j  ,k2))+              &
     &                           tl_cff4*(cff4*dTdz(i,j,k2)-            &
     &                                    dTde(i,j+1,k1))+              &
     &                           cff1*(tl_cff1*dTdz(i,j,k2)+            &
     &                                 cff1*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j  ,k1))+              &
     &                           cff2*(tl_cff2*dTdz(i,j,k2)+            &
     &                                 cff2*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j+1,k2))+              &
     &                           cff3*(tl_cff3*dTdz(i,j,k2)+            &
     &                                 cff3*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j  ,k2))+              &
     &                           cff4*(tl_cff4*dTdz(i,j,k2)+            &
     &                                 cff4*tl_dTdz(i,j,k2)-            &
     &                                 tl_dTde(i,j+1,k1)))
                END DO
              END DO
            END IF
!
!  Time-step biharmonic, geopotential diffusion term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
!>              cff=dt(ng)*pm(i,j)*pn(i,j)*                             &
!>   &                     (FX(i+1,j  )-FX(i,j)+                        &
!>   &                      FE(i  ,j+1)-FE(i,j))+                       &
!>   &              dt(ng)*(FS(i,j,k2)-FS(i,j,k1))
!>
                tl_cff=dt(ng)*pm(i,j)*pn(i,j)*                          &
     &                        (tl_FX(i+1,j)-tl_FX(i,j)+                 &
     &                         tl_FE(i,j+1)-tl_FE(i,j))+                &
     &                 dt(ng)*(tl_FS(i,j,k2)-tl_FS(i,j,k1))
!>              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff
!>
                tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff
#ifdef DIAGNOSTICS_TS
!!              DiaTwrk(i,j,k,itrc,iThdif)=-cff
#endif
              END DO
            END DO
          END IF
        END DO K_LOOP2
      END DO T_LOOP

      RETURN
      END SUBROUTINE tl_t3dmix4_tile
