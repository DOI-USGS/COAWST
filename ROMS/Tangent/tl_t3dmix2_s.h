      SUBROUTINE tl_t3dmix2 (ng, tile)
!
!svn $Id: tl_t3dmix2_s.h 751 2015-01-07 22:56:36Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2015 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes tangent linear horizontal harmonic mixing  !
!  of tracers along S-coordinate levels surfaces.                      !
!                                                                      !
!  BASIC STATE variables needed:  t, Hz                                !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef CLIMA_TS_MIX
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
      CALL wclock_on (ng, iTLM, 24)
#endif
      CALL tl_t3dmix2_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % pmon_u,                          &
     &                      GRID(ng) % pnom_v,                          &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      MIXING(ng) % diff2,                         &
#ifdef CLIMA_TS_MIX
     &                      CLIMA(ng) % tclm,                           &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                      DIAGS(ng) % DiaTwrk,                        &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)
#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 24)
#endif
      RETURN
      END SUBROUTINE tl_t3dmix2
!
!***********************************************************************
      SUBROUTINE tl_t3dmix2_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nnew,                           &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
     &                            Hz, tl_Hz,                            &
     &                            pmon_u, pnom_v, pm, pn,               &
     &                            diff2,                                &
#ifdef CLIMA_TS_MIX
     &                            tclm,                                 &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                            DiaTwrk,                              &
#endif
     &                            t, tl_t)
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
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
# ifdef CLIMA_TS_MIX
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
# ifdef CLIMA_TS_MIX
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: cff, cff1, tl_cff, tl_cff1

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute tangent linear horizontal harmonic diffusion along constant
!  S-surfaces.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*          &
     &            pmon_u(i,j)
!>            FX(i,j)=cff*
!>   &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
#if defined CLIMA_TS_MIX
!>   &                ((t(i  ,j,k,nrhs,itrc)-tclm(i  ,j,k,itrc))-       &
!>   &                 (t(i-1,j,k,nrhs,itrc)-tclm(i-1,j,k,itrc)))
#else
!>   &                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
#endif
!>
              tl_FX(i,j)=cff*                                           &
#if defined CLIMA_TS_MIX
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    ((t(i  ,j,k,nrhs,itrc)-tclm(i  ,j,k,itrc))-   &
     &                     (t(i-1,j,k,nrhs,itrc)-tclm(i-1,j,k,itrc)))+  &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_t(i  ,j,k,nrhs,itrc)-                     &
     &                     tl_t(i-1,j,k,nrhs,itrc)))
#else
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (t(i  ,j,k,nrhs,itrc)-                        &
     &                     t(i-1,j,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_t(i  ,j,k,nrhs,itrc)-                     &
     &                     tl_t(i-1,j,k,nrhs,itrc)))
#endif
#ifdef MASKING
!>            FX(i,j)=FX(i,j)*umask(i,j)
!>
              tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
#endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*          &
     &            pnom_v(i,j)
!>            FE(i,j)=cff*                                              &
!>   &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
#if defined CLIMA_TS_MIX
!>   &                ((t(i,j  ,k,nrhs,itrc)-tclm(i,j  ,k,itrc))-       &
!>   &                 (t(i,j-1,k,nrhs,itrc)-tclm(i,j-1,k,itrc)))
#else
!>   &                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
#endif
!>
              tl_FE(i,j)=cff*                                           &
#if defined CLIMA_TS_MIX
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    ((t(i,j  ,k,nrhs,itrc)-tclm(i,j  ,k,itrc))-   &
     &                     (t(i,j-1,k,nrhs,itrc)-tclm(i,j-1,k,itrc)))+
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_t(i,j  ,k,nrhs,itrc)-                     &
     &                     tl_t(i,j-1,k,nrhs,itrc)))
#else
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (t(i,j  ,k,nrhs,itrc)-                        &
     &                     t(i,j-1,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_t(i,j  ,k,nrhs,itrc)-                     &
     &                     tl_t(i,j-1,k,nrhs,itrc)))
#endif
#ifdef MASKING
!>            FE(i,j)=FE(i,j)*vmask(i,j)
!>
              tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
#endif
            END DO
          END DO
!
! Time-step harmonic, S-surfaces diffusion term (m Tunits).
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
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
#ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE tl_t3dmix2_tile
