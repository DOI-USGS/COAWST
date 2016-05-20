      SUBROUTINE t3dmix4 (ng, tile)
!
!svn $Id: t3dmix4_s.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine computes horizontal biharmonic mixing of tracers    !
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
      CALL wclock_on (ng, iNLM, 27)
#endif
      CALL t3dmix4_tile (ng, tile,                                      &
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
# ifdef TS_U3ADV_SPLIT
     &                   MIXING(ng) % diff3d_u,                         &
     &                   MIXING(ng) % diff3d_v,                         &
# else
     &                   MIXING(ng) % diff3d_r,                         &
# endif
#else
     &                   MIXING(ng) % diff4,                            &
#endif
#ifdef TS_MIX_CLIMA
     &                   CLIMA(ng) % tclm,                              &
#endif
#ifdef DIAGNOSTICS_TS
     &                   DIAGS(ng) % DiaTwrk,                           &
#endif
     &                   OCEAN(ng) % t)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 27)
#endif

      RETURN
      END SUBROUTINE t3dmix4
!
!***********************************************************************
      SUBROUTINE t3dmix4_tile (ng, tile,                                &
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
# ifdef TS_U3ADV_SPLIT
     &                         diff3d_u, diff3d_v,                      &
# else
     &                         diff3d_r,                                &
# endif
#else
     &                         diff4,                                   &
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
      USE mod_ncparam
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
#  ifdef TS_U3ADV_SPLIT
      real(r8), intent(in) :: diff3d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: diff3d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: diff3d_r(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:,LBj:,:)
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
#  ifdef TS_U3ADV_SPLIT
      real(r8), intent(in) :: diff3d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: diff3d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
# else
      real(r8), intent(in) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
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
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ibt, itrc, j, k

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapT

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
#ifdef OFFLINE_BIOLOGY
      DO ibt=1,NBT
        itrc=idbio(ibt)
#else
      DO itrc=1,NT(ng)
#endif
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax+1
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
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
#ifdef WET_DRY
              cff=cff*umask_diff(i,j)
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
     &                  (t(i  ,j,k,nrhs,itrc)-                          &
     &                   t(i-1,j,k,nrhs,itrc))
              END IF
#else
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (t(i  ,j,k,nrhs,itrc)-                            &
     &                 t(i-1,j,k,nrhs,itrc))
#endif
            END DO
          END DO
          DO j=Jmin,Jmax+1
            DO i=Imin,Imax
#ifdef DIFF_3DCOEF
# ifdef TS_U3ADV_SPLIT
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
#ifdef WET_DRY
              cff=cff*vmask_diff(i,j)
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
     &                  (t(i,j  ,k,nrhs,itrc)-                          &
     &                   t(i,j-1,k,nrhs,itrc))
              END IF
#else
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (t(i,j  ,k,nrhs,itrc)-                            &
     &                 t(i,j-1,k,nrhs,itrc))
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
              LapT(i,j)=pm(i,j)*pn(i,j)*cff*                            &
     &                  (FX(i+1,j)-FX(i,j)+                             &
     &                   FE(i,j+1)-FE(i,j))
            END DO
          END DO
!
!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first harmonic operator.
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              IF (LBC(iwest,isTvar(itrc),ng)%closed) THEN
                DO j=Jmin,Jmax
                  LapT(Istr-1,j)=0.0_r8
                END DO
              ELSE
                DO j=Jmin,Jmax
                  LapT(Istr-1,j)=LapT(Istr,j)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              IF (LBC(ieast,isTvar(itrc),ng)%closed) THEN
                DO j=Jmin,Jmax
                  LapT(Iend+1,j)=0.0_r8
                END DO
              ELSE
                DO j=Jmin,Jmax
                  LapT(Iend+1,j)=LapT(Iend,j)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              IF (LBC(isouth,isTvar(itrc),ng)%closed) THEN
                DO i=Imin,Imax
                  LapT(i,Jstr-1)=0.0_r8
                END DO
              ELSE
                DO i=Imin,Imax
                  LapT(i,Jstr-1)=LapT(i,Jstr)
                END DO
              END IF
            END IF
          END IF
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              IF (LBC(inorth,isTvar(itrc),ng)%closed) THEN
                DO i=Imin,Imax
                  LapT(i,Jend+1)=0.0_r8
                END DO
              ELSE
                DO i=Imin,Imax
                  LapT(i,Jend+1)=LapT(i,Jend)
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
# ifdef TS_U3ADV_SPLIT
              cff=0.5_r8*diff3d_u(i,j,k)*pmon_u(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*          &
     &            pmon_u(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#endif
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (LapT(i,j)-LapT(i-1,j))
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
# ifdef TS_U3ADV_SPLIT
              cff=0.5_r8*diff3d_v(i,j,k)*pnom_v(i,j)
# else
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*          &
     &            pnom_v(i,j)
# endif
#else
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*          &
     &            pnom_v(i,j)
#endif
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (LapT(i,j)-LapT(i,j-1))
#ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#endif
#ifdef WET_DRY
              FE(i,j)=FE(i,j)*vmask_diff(i,j)
#endif
            END DO
          END DO
!
!  Time-step biharmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j  )-FX(i,j))
              cff2=cff*(FE(i  ,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
#ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxdif)=-cff1
              DiaTwrk(i,j,k,itrc,iTydif)=-cff2
              DiaTwrk(i,j,k,itrc,iThdif)=-cff3
#endif
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE t3dmix4_tile
