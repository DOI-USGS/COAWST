#ifdef EW_PERIODIC
# define I_RANGE Istr-1,Iend+1
#else
# define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#endif
#ifdef NS_PERIODIC
# define J_RANGE Jstr-1,Jend+1
#else
# define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#endif
#define MIX_STABILITY

      SUBROUTINE t3dmix4 (ng, tile)
!
!svn $Id: t3dmix4_s.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
     &                   nrhs(ng), nstp(ng), nnew(ng),                  &
#ifdef MASKING
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
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
     &                         nrhs, nstp, nnew,                        &
#ifdef MASKING
     &                         umask, vmask,                            &
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
#ifdef DIAGNOSTICS_TS
     &                         DiaTwrk,                                 &
#endif
     &                         t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
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
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
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
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: cff, cff1

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: LapT

#include "set_bounds.h"
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic diffusion along constant S-surfaces.
!  The biharmonic operator is computed by applying the harmonic
!  operator twice.
#ifdef MIX_STABILITY
!  In order to increase stability, the biharmonic operator is applied
!  as: 3/4 t(:,:,:,nrhs,:) + 1/4 t(:,:,:,nstp,:).
#endif
!-----------------------------------------------------------------------
!
!  Compute horizontal tracer flux in the XI- and ETA-directions.
!
          DO j=J_RANGE
            DO i=I_RANGE+1
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
              FX(i,j)=cff*(Hz(i,j,k)+Hz(i-1,j,k))*                      &
#ifdef MIX_STABILITY
     &                (0.75_r8*(t(i  ,j,k,nrhs,itrc)-                   &
     &                          t(i-1,j,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i  ,j,k,nstp,itrc)-                   &
     &                          t(i-1,j,k,nstp,itrc)))
#else
     &                (t(i  ,j,k,nrhs,itrc)-                            &
     &                 t(i-1,j,k,nrhs,itrc))
#endif
            END DO
          END DO
          DO j=J_RANGE+1
            DO i=I_RANGE
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
              FE(i,j)=cff*(Hz(i,j,k)+Hz(i,j-1,k))*                      &
#ifdef MIX_STABILITY
     &                (0.75_r8*(t(i,j  ,k,nrhs,itrc)-                   &
     &                          t(i,j-1,k,nrhs,itrc))+                  &
     &                 0.25_r8*(t(i,j  ,k,nstp,itrc)-                   &
     &                          t(i,j-1,k,nstp,itrc)))
#else
     &                (t(i,j  ,k,nrhs,itrc)-                            &
     &                 t(i,j-1,k,nrhs,itrc))
#endif
            END DO
          END DO
!
!  Compute first harmonic operator and multiply by the metrics of the
!  second harmonic operator.
!
          DO j=J_RANGE
            DO i=I_RANGE
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
#ifndef EW_PERIODIC
          IF (WESTERN_EDGE) THEN
            DO j=J_RANGE
# ifdef WESTERN_WALL
              LapT(Istr-1,j)=0.0_r8
# else
              LapT(Istr-1,j)=LapT(Istr,j)
# endif
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=J_RANGE
# ifdef EASTERN_WALL
              LapT(Iend+1,j)=0.0_r8
# else
              LapT(Iend+1,j)=LapT(Iend,j)
# endif
            END DO
          END IF
#endif
#ifndef NS_PERIODIC
          IF (SOUTHERN_EDGE) THEN
            DO i=I_RANGE
# ifdef SOUTHERN_WALL
              LapT(i,Jstr-1)=0.0_r8
# else
              LapT(i,Jstr-1)=LapT(i,Jstr)
# endif
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO i=I_RANGE
# ifdef NORTHERN_WALL
              LapT(i,Jend+1)=0.0_r8
# else
              LapT(i,Jend+1)=LapT(i,Jend)
# endif
            END DO
          END IF
#endif
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
            END DO
          END DO
!
!  Time-step biharmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)*                               &
     &                   (FX(i+1,j)-FX(i,j)+                            &
     &                    FE(i,j+1)-FE(i,j))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff
#ifdef TS_MPDATA
              cff1=1.0_r8/Hz(i,j,k)
              t(i,j,k,3,itrc)=cff1*t(i,j,k,nnew,itrc)
#endif
#ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iThdif)=-cff
#endif
            END DO
          END DO
        END DO
      END DO
#undef I_RANGE
#undef J_RANGE
      RETURN
      END SUBROUTINE t3dmix4_tile
