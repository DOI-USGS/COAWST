      SUBROUTINE rp_t3dmix4 (ng, tile)
!
!svn $Id: rp_t3dmix4_s.h 995 2020-01-10 04:01:28Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2020 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes representers tangent linear horizontal     !
!  biharmonic mixing of tracers along S-coordinate levels surfaces.    !
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
      CALL wclock_on (ng, iRPM, 27, __LINE__, __FILE__)
#endif
      CALL rp_t3dmix4_tile (ng, tile,                                   &
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
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % pmon_u,                          &
     &                      GRID(ng) % pnom_v,                          &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
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
      CALL wclock_off (ng, iRPM, 27, __LINE__, __FILE__)
#endif

      RETURN
      END SUBROUTINE rp_t3dmix4
!
!***********************************************************************
      SUBROUTINE rp_t3dmix4_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
#ifdef WET_DRY_NOT_YET
     &                            umask_wet, vmask_wet,                 &
#endif
     &                            Hz, tl_Hz,                            &
     &                            pmon_u, pnom_v, pm, pn,               &
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
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
#  ifdef TS_U3ADV_SPLIT_NOT_YET
      real(r8), intent(in) :: diff3d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: diff3d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: diff3d_r(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:,LBj:,:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
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
#  ifdef TS_U3ADV_SPLIT_NOT_YET
      real(r8), intent(in) :: diff3d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: diff3d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
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
      integer :: i, itrc, j, k

      real(r8) :: cff, cff1, tl_cff, tl_cff1

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapT

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LapT

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic diffusion along constant S-surfaces.
!  The biharmonic operator is computed by applying the harmonic
!  operator twice.
#ifdef TS_MIX_STABILITY
!  In order to increase stability, the biharmonic operator is applied
!  as: 3/4 t(:,:,:,nrhs,:) + 1/4 t(:,:,:,nstp,:).
#endif
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
!  Compute horizontal tracer flux in the XI- and ETA-directions.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax+1
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
              cff=0.5_r8*diff3d_u(i,j,k)*pmon_u(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*          &
     &            pmon_u(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#endif
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              cff=cff*umask_wet(i,j)
#endif
#if defined TS_MIX_STABILITY
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (0.75_r8*(t(i  ,j,k,nrhs,itrc)-                   &
     &                          t(i-1,j,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i  ,j,k,nstp,itrc)-                   &
     &                          t(i-1,j,k,nstp,itrc)))
              tl_FX(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (0.75_r8*(t(i  ,j,k,nrhs,itrc)-               &
     &                              t(i-1,j,k,nrhs,itrc))+              &
     &                     0.25_r8*(t(i  ,j,k,nstp,itrc)-               &
     &                              t(i-1,j,k,nstp,itrc)))+             &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (0.75_r8*(tl_t(i  ,j,k,nrhs,itrc)-            &
     &                              tl_t(i-1,j,k,nrhs,itrc))+           &
     &                     0.25_r8*(tl_t(i  ,j,k,nstp,itrc)-            &
     &                              tl_t(i-1,j,k,nstp,itrc))))-         &
# ifdef TL_IOMS
     &                   FX(i,j)
# endif
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  ((t(i  ,j,k,nrhs,itrc)-tclm(i  ,j,k,itrc))-     &
     &                   (t(i-1,j,k,nrhs,itrc)-tclm(i-1,j,k,itrc)))
                tl_FX(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      ((t(i  ,j,k,nrhs,itrc)-                     &
     &                        tclm(i  ,j,k,itrc))-                      &
     &                       (t(i-1,j,k,nrhs,itrc)-                     &
     &                        tclm(i-1,j,k,itrc)))+                     &
     &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      (tl_t(i  ,j,k,nrhs,itrc)-                   &
     &                       tl_t(i-1,j,k,nrhs,itrc)))-                 &
# ifdef TL_IOMS
     &                     FX(i,j)
# endif
              ELSE
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (t(i  ,j,k,nrhs,itrc)-                          &
     &                   t(i-1,j,k,nrhs,itrc))
                tl_FX(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      (t(i  ,j,k,nrhs,itrc)-                      &
     &                       t(i-1,j,k,nrhs,itrc))+                     &
     &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      (tl_t(i  ,j,k,nrhs,itrc)-                   &
     &                       tl_t(i-1,j,k,nrhs,itrc)))-                 &
# ifdef TL_IOMS
     &                     FX(i,j)
# endif
              END IF
#else
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (t(i  ,j,k,nrhs,itrc)-                            &
     &                 t(i-1,j,k,nrhs,itrc))
              tl_FX(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (t(i  ,j,k,nrhs,itrc)-                        &
     &                     t(i-1,j,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_t(i  ,j,k,nrhs,itrc)-                     &
     &                     tl_t(i-1,j,k,nrhs,itrc)))-                   &
# ifdef TL_IOMS
     &                   FX(i,j)
# endif
#endif
            END DO
          END DO
          DO j=Jmin,Jmax+1
            DO i=Imin,Imax
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
              cff=0.5_r8*diff3d_v(i,j,k)*pnom_v(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*          &
     &            pnom_v(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*          &
     &            pnom_v(i,j)
#endif
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              cff=cff*vmask_wet(i,j)
#endif
#if defined TS_MIX_STABILITY
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (0.75_r8*(t(i,j  ,k,nrhs,itrc)-                   &
     &                          t(i,j-1,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i,j  ,k,nstp,itrc)-                   &
     &                          t(i,j-1,k,nstp,itrc)))
              tl_FE(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (0.75_r8*(t(i,j  ,k,nrhs,itrc)-               &
     &                              t(i,j-1,k,nrhs,itrc))+              &
     &                     0.25_r8*(t(i,j  ,k,nstp,itrc)-               &
     &                              t(i,j-1,k,nstp,itrc)))+             &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (0.75_r8*(tl_t(i,j  ,k,nrhs,itrc)-            &
     &                              tl_t(i,j-1,k,nrhs,itrc))+           &
     &                     0.25_r8*(tl_t(i,j  ,k,nstp,itrc)-            &
     &                              tl_t(i,j-1,k,nstp,itrc))))-         &
# ifdef TL_IOMS
     &                   FE(i,j)
# endif
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  ((t(i,j  ,k,nrhs,itrc)-tclm(i,j  ,k,itrc))-     &
     &                   (t(i,j-1,k,nrhs,itrc)-tclm(i,j-1,k,itrc)))
                tl_FE(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      ((t(i,j  ,k,nrhs,itrc)-                     &
     &                        tclm(i,j  ,k,itrc))-                      &
     &                       (t(i,j-1,k,nrhs,itrc)-                     &
     &                        tclm(i,j-1,k,itrc)))+                     &
     &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      (tl_t(i,j  ,k,nrhs,itrc)-                   &
     &                       tl_t(i,j-1,k,nrhs,itrc)))-                 &
# ifdef TL_IOMS
     &                     FE(i,j)
# endif
              ELSE
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (t(i,j  ,k,nrhs,itrc)-                          &
     &                   t(i,j-1,k,nrhs,itrc))
                tl_FE(i,j)=cff*                                         &
     &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      (t(i,j  ,k,nrhs,itrc)-                      &
     &                       t(i,j-1,k,nrhs,itrc))+                     &
     &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      (tl_t(i,j  ,k,nrhs,itrc)-                   &
     &                       tl_t(i,j-1,k,nrhs,itrc)))-                 &
# ifdef TL_IOMS
     &                     FE(i,j)
# endif
              END IF
#else
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (t(i,j  ,k,nrhs,itrc)-                            &
     &                 t(i,j-1,k,nrhs,itrc))
              tl_FE(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (t(i,j  ,k,nrhs,itrc)-                        &
     &                     t(i,j-1,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_t(i,j  ,k,nrhs,itrc)-                     &
     &                     tl_t(i,j-1,k,nrhs,itrc)))-                   &
# ifdef TL_IOMS
     &                   FE(i,j)
# endif
#endif
            END DO
          END DO
!
!  Compute first harmonic operator and multiply by the metrics of the
!  second harmonic operator.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              cff=1.0_r8/Hz(i,j,k)
              tl_cff=-cff*cff*tl_Hz(i,j,k)+                             &
#ifdef TL_IOMS
     &               2.0_r8*cff
#endif
              LapT(i,j)=pm(i,j)*pn(i,j)*cff*                            &
     &                  (FX(i+1,j)-FX(i,j)+                             &
     &                   FE(i,j+1)-FE(i,j))
              tl_LapT(i,j)=pm(i,j)*pn(i,j)*                             &
     &                     (tl_cff*                                     &
     &                      (FX(i+1,j)-FX(i,j)+                         &
     &                       FE(i,j+1)-FE(i,j))+                        &
     &                      cff*                                        &
     &                      (tl_FX(i+1,j)-tl_FX(i,j)+                   &
     &                       tl_FE(i,j+1)-tl_FE(i,j)))-                 &
#ifdef TL_IOMS
     &                     LapT(i,j)
#endif
            END DO
          END DO
!
!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first harmonic operator.
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              IF (tl_LBC(iwest,isTvar(itrc),ng)%closed) THEN
                DO j=Jmin,Jmax
                  LapT(Istr-1,j)=0.0_r8
                  tl_LapT(Istr-1,j)=0.0_r8
                END DO
              ELSE
                DO j=Jmin,Jmax
                  LapT(Istr-1,j)=LapT(Istr,j)
                  tl_LapT(Istr-1,j)=tl_LapT(Istr,j)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              IF (tl_LBC(ieast,isTvar(itrc),ng)%closed) THEN
                DO j=Jmin,Jmax
                  LapT(Iend+1,j)=0.0_r8
                  tl_LapT(Iend+1,j)=0.0_r8
                END DO
              ELSE
                DO j=Jmin,Jmax
                  LapT(Iend+1,j)=LapT(Iend,j)
                  tl_LapT(Iend+1,j)=tl_LapT(Iend,j)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              IF (tl_LBC(isouth,isTvar(itrc),ng)%closed) THEN
                DO i=Imin,Imax
                  LapT(i,Jstr-1)=0.0_r8
                  tl_LapT(i,Jstr-1)=0.0_r8
                END DO
              ELSE
                DO i=Imin,Imax
                  LapT(i,Jstr-1)=LapT(i,Jstr)
                  tl_LapT(i,Jstr-1)=tl_LapT(i,Jstr)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              IF (tl_LBC(inorth,isTvar(itrc),ng)%closed) THEN
                DO i=Imin,Imax
                  LapT(i,Jend+1)=0.0_r8
                  tl_LapT(i,Jend+1)=0.0_r8
                END DO
              ELSE
                DO i=Imin,Imax
                  LapT(i,Jend+1)=LapT(i,Jend)
                  tl_LapT(i,Jend+1)=tl_LapT(i,Jend)
                END DO
              END IF
            END IF
          END IF
!
!  Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
              cff=0.5_r8*diff3d_u(i,j,k)*pmon_u(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*          &
     &            pmon_u(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#endif
!>            FX(i,j)=cff*                                              &
!>   &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
!>   &                (LapT(i,j)-LapT(i-1,j))
!>
              tl_FX(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (LapT(i,j)-LapT(i-1,j))+                      &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_LapT(i,j)-tl_LapT(i-1,j)))-               &
# ifdef TL_IOMS
     &                   cff*                                           &
     &                   (Hz(i,j,k)+Hz(i-1,j,k))*                       &
     &                   (LapT(i,j)-LapT(i-1,j))
# endif
#ifdef MASKING
!>            FX(i,j)=FX(i,j)*umask(i,j)
!>
              tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              FX(i,j)=FX(i,j)*umask_wet(i,j)
#endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT_NOT_YET
              cff=0.5_r8*diff3d_v(i,j,k)*pnom_v(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*          &
     &            pnom_v(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*          &
     &            pnom_v(i,j)
#endif
!>            FE(i,j)=cff*                                              &
!>   &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
!>   &                (LapT(i,j)-LapT(i,j-1))
!>
              tl_FE(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (LapT(i,j)-LapT(i,j-1))+                      &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_LapT(i,j)-tl_LapT(i,j-1)))-               &
#ifdef TL_IOMS
     &                   cff*                                           &
     &                   (Hz(i,j,k)+Hz(i,j-1,k))*                       &
     &                   (LapT(i,j)-LapT(i,j-1))
#endif
#ifdef MASKING
!>            FE(i,j)=FE(i,j)*vmask(i,j)
!>
              tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              FE(i,j)=FE(i,j)*vmask_wet(i,j)
#endif
            END DO
          END DO
!
!  Time-step biharmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            cff=dt(ng)*pm(i,j)*pn(i,j)*                               &
!>   &                   (FX(i+1,j)-FX(i,j)+                            &
!>   &                    FE(i,j+1)-FE(i,j))
!>
              tl_cff=dt(ng)*pm(i,j)*pn(i,j)*                            &
     &                      (tl_FX(i+1,j)-tl_FX(i,j)+                   &
     &                       tl_FE(i,j+1)-tl_FE(i,j))
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff
#ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iThdif)=-cff
#endif
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE rp_t3dmix4_tile
