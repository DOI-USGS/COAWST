      SUBROUTINE t3dmix2 (ng, tile)
!
!svn $Id: t3dmix2_s.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine computes horizontal harmonic mixing of tracers      !
!  along S-coordinate levels surfaces.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef TS_MIX_CLIMA
      USE mod_clima
#endif
#ifdef DIAGNOSTICS_TS
      USE mod_diags
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
      CALL wclock_on (ng, iNLM, 24)
#endif
      CALL t3dmix2_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng), nstp(ng), nnew(ng),                  &
#ifdef MASKING
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#endif
#ifdef WET_DRY
     &                   GRID(ng) % umask_diff,                         &
     &                   GRID(ng) % vmask_diff,                         &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % pmon_u,                             &
     &                   GRID(ng) % pnom_v,                             &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
#ifdef DIFF_3DCOEF
     &                   MIXING(ng) % diff3d_r,                         &
#else
     &                   MIXING(ng) % diff2,                            &
#endif
#ifdef TS_MIX_CLIMA
     &                   CLIMA(ng) % tclm,                              &
#endif
#ifdef DIAGNOSTICS_TS
     &                   DIAGS(ng) % DiaTwrk,                           &
#endif
     &                   OCEAN(ng) % t)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 24)
#endif
      RETURN
      END SUBROUTINE t3dmix2
!
!***********************************************************************
      SUBROUTINE t3dmix2_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs, nstp, nnew,                        &
#ifdef MASKING
     &                         umask, vmask,                            &
#endif
#ifdef WET_DRY
     &                         umask_diff, vmask_diff,                  &
#endif
     &                         Hz, pmon_u, pnom_v, pm, pn,              &
#ifdef DIFF_3DCOEF
     &                         diff3d_r,                                &
#else
     &                         diff2,                                   &
#endif
#ifdef TS_MIX_CLIMA
     &                         tclm,                                    &
#endif
#ifdef DIAGNOSTICS_TS
     &                         DiaTwrk,                                 &
#endif
     &                         t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
#ifdef OFFLINE_BIOLOGY
      USE mod_biology
#endif
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
# ifdef WET_DRY
      real(r8), intent(in) :: umask_diff(LBi:,LBj:)
      real(r8), intent(in) :: vmask_diff(LBi:,LBj:)
# endif
# ifdef DIFF_3DCOEF
      real(r8), intent(in) :: diff3d_r(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: umask_diff(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_diff(LBi:UBi,LBj:UBj)
# endif
# ifdef DIFF_3DCOEF
      real(r8), intent(in) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, ibt, itrc, j, k

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic diffusion along constant S-surfaces.
#ifdef TS_MIX_STABILITY
!  In order to increase stability, the harmonic operator is applied
!  as: 3/4 t(:,:,:,nrhs,:) + 1/4 t(:,:,:,nstp,:).
#endif
!-----------------------------------------------------------------------
!
#ifdef OFFLINE_BIOLOGY
      DO ibt=1,NBT
        itrc=idbio(ibt)
#else
      DO itrc=1,NT(ng)
#endif
        DO k=1,N(ng)
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#ifdef DIFF_3DCOEF
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*          &
     &            pmon_u(i,j)
#else
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#endif
#if defined TS_MIX_STABILITY
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (0.75_r8*(t(i  ,j,k,nrhs,itrc)-                   &
     &                          t(i-1,j,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i  ,j,k,nstp,itrc)-                   &
     &                          t(i-1,j,k,nstp,itrc)))
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  ((t(i  ,j,k,nrhs,itrc)-tclm(i  ,j,k,itrc))-     &
     &                   (t(i-1,j,k,nrhs,itrc)-tclm(i-1,j,k,itrc)))
              ELSE
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
              END IF
#else
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
#endif
#ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#endif
#ifdef WET_DRY
              FX(i,j)=FX(i,j)*umask_diff(i,j)
#endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#ifdef DIFF_3DCOEF
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*          &
     &            pnom_v(i,j)
#else
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*          &
     &            pnom_v(i,j)
#endif
#if defined TS_MIX_STABILITY
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (0.75_r8*(t(i,j  ,k,nrhs,itrc)-                   &
     &                          t(i,j-1,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i,j  ,k,nstp,itrc)-                   &
     &                          t(i,j-1,k,nstp,itrc)))
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  ((t(i,j  ,k,nrhs,itrc)-tclm(i,j  ,k,itrc))-     &
     &                   (t(i,j-1,k,nrhs,itrc)-tclm(i,j-1,k,itrc)))
              ELSE
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
              END IF
#else
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
#endif
#ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#endif
#ifdef WET_DRY
              FE(i,j)=FE(i,j)*vmask_diff(i,j)
#endif
            END DO
          END DO
!
! Time-step harmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j  )-FX(i,j))
              cff2=cff*(FE(i  ,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff3
#ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxdif)=cff1
              DiaTwrk(i,j,k,itrc,iTydif)=cff2
              DiaTwrk(i,j,k,itrc,iThdif)=cff3
#endif
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE t3dmix2_tile
