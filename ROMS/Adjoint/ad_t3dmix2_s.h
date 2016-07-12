      SUBROUTINE ad_t3dmix2 (ng, tile)
!
!svn $Id: ad_t3dmix2_s.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes adjoint horizontal harmonic mixing of      !
!  tracers along S-coordinate levels surfaces.                         !
!                                                                      !
!  BASIC STATE variables needed:  t, Hz                                !
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
      CALL wclock_on (ng, iADM, 24)
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
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
     &                      GRID(ng) % pmon_u,                          &
     &                      GRID(ng) % pnom_v,                          &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
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
      CALL wclock_off (ng, iADM, 24)
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
     &                            Hz, ad_Hz,                            &
     &                            pmon_u, pnom_v, pm, pn,               &
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
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
# ifdef TS_MIX_CLIMA
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
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
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: cff, cff1
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4, ad_cff, ad_cff1

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8

      ad_FE(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS,JminS:JmaxS)=0.0_r8
!
!----------------------------------------------------------------------
!  Compute adjoint horizontal harmonic diffusion along constant
!  S-surfaces.
!----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
! Time-step harmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
#ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
!>            tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
!>
              ad_cff=ad_cff+ad_t(i,j,k,nnew,itrc)
!>            tl_cff=dt(ng)*pm(i,j)*pn(i,j)*                            &
!>   &                      (tl_FX(i+1,j)-tl_FX(i,j)+                   &
!>   &                       tl_FE(i,j+1)-tl_FE(i,j))
!>
              adfac=dt(ng)*pm(i,j)*pn(i,j)*ad_cff
              ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
              ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
              ad_FE(i,j  )=ad_FE(i,j  )-adfac
              ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#ifdef DIFF_3DCOEF
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))*          &
     &            pnom_v(i,j)
#else
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*          &
     &            pnom_v(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              FE(i,j)=FE(i,j)*vmask_wet(i,j)
#endif
#ifdef MASKING
!>            tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
!>
              ad_FE(i,j)=ad_FE(i,j)*vmask(i,j)
#endif
#if defined TS_MIX_STABILITY
!>            tl_FE(i,j)=cff*                                           &
!>   &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
!>   &                    (0.75_r8*(t(i,j  ,k,nrhs,itrc)-               &
!>   &                              t(i,j-1,k,nrhs,itrc))+              &
!>   &                     0.25_r8*(t(i,j  ,k,nstp,itrc)-               &
!>   &                              t(i,j-1,k,nstp,itrc)))+             &
!>   &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
!>   &                    (0.75_r8*(tl_t(i,j  ,k,nrhs,itrc)-            &
!>   &                              tl_t(i,j-1,k,nrhs,itrc))+           &
!>   &                     0.25_r8*(tl_t(i,j  ,k,nstp,itrc)-            &
!>   &                              tl_t(i,j-1,k,nstp,itrc))))
!>
              adfac=cff*ad_FE(i,j)
              adfac1=adfac*(0.75_r8*(t(i,j  ,k,nrhs,itrc)-              &
     &                               t(i,j-1,k,nrhs,itrc))+             &
     &                      0.25_r8*(t(i,j  ,k,nstp,itrc)-              &
     &                               t(i,j-1,k,nstp,itrc)))
              adfac2=adfac*(Hz(i,j,k)+Hz(i,j-1,k))
              adfac3=adfac2*0.75_r8
              adfac4=adfac2*0.25_r8
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
              ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-adfac3
              ad_t(i,j  ,k,nrhs,itrc)=ad_t(i,j  ,k,nrhs,itrc)+adfac3
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)-adfac4
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+adfac4
              ad_FE(i,j)=0.0_r8
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
!>              tl_FE(i,j)=cff*                                         &
!>   &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
!>   &                      ((t(i,j  ,k,nrhs,itrc)-                     &
!>   &                        tclm(i,j  ,k,itrc))-                      &
!>   &                       (t(i,j-1,k,nrhs,itrc)-                     &
!>   &                        tclm(i,j-1,k,itrc)))+                     &
!>   &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
!>   &                      (tl_t(i,j  ,k,nrhs,itrc)-                   &
!>   &                       tl_t(i,j-1,k,nrhs,itrc)))
!>
                adfac=cff*ad_FE(i,j)
                adfac1=adfac*((t(i,j  ,k,nrhs,itrc)-                    &
     &                         tclm(i,j  ,k,itrc))-                     &
     &                        (t(i,j-1,k,nrhs,itrc)-                    &
     &                         tclm(i,j-1,k,itrc)))
                adfac2=adfac*(Hz(i,j,k)+Hz(i,j-1,k))
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-adfac2
                ad_t(i,j  ,k,nrhs,itrc)=ad_t(i,j  ,k,nrhs,itrc)+adfac2
                ad_FE(i,j)=0.0_r8
              ELSE
!>              tl_FE(i,j)=cff*                                         &
!>   &                     ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
!>   &                      (t(i,j  ,k,nrhs,itrc)-                      &
!>   &                       t(i,j-1,k,nrhs,itrc))+                     &
!>   &                      (Hz(i,j,k)+Hz(i,j-1,k))*                    &
!>   &                      (tl_t(i,j  ,k,nrhs,itrc)-                   &
!>   &                       tl_t(i,j-1,k,nrhs,itrc)))
!>
                adfac=cff*ad_FE(i,j)
                adfac1=adfac*(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
                adfac2=adfac*(Hz(i,j,k)+Hz(i,j-1,k))
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-adfac2
                ad_t(i,j  ,k,nrhs,itrc)=ad_t(i,j  ,k,nrhs,itrc)+adfac2
                ad_FE(i,j)=0.0_r8
              END IF
#else
!>            tl_FE(i,j)=cff*                                           &
!>   &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
!>   &                    (t(i,j  ,k,nrhs,itrc)-                        &
!>   &                     t(i,j-1,k,nrhs,itrc))+                       &
!>   &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
!>   &                    (tl_t(i,j  ,k,nrhs,itrc)-                     &
!>   &                     tl_t(i,j-1,k,nrhs,itrc)))
!>
              adfac=cff*ad_FE(i,j)
              adfac1=adfac*(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
              adfac2=adfac*(Hz(i,j,k)+Hz(i,j-1,k))
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
              ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-adfac2
              ad_t(i,j  ,k,nrhs,itrc)=ad_t(i,j  ,k,nrhs,itrc)+adfac2
              ad_FE(i,j)=0.0_r8
#endif
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#ifdef DIFF_3DCOEF
              cff=0.25_r8*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))*          &
     &            pmon_u(i,j)
#else
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#endif
#ifdef WET_DRY_NOT_YET
              FX(i,j)=FX(i,j)*umask_wet(i,j)
#endif
#ifdef MASKING
!>            tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
!>
              ad_FX(i,j)=ad_FX(i,j)*umask(i,j)
#endif
#if defined TS_MIX_STABILITY
!>            tl_FX(i,j)=cff*                                           &
!>   &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
!>   &                    (0.75_r8*(t(i  ,j,k,nrhs,itrc)-               &
!>   &                              t(i-1,j,k,nrhs,itrc))+              &
!>   &                     0.25_r8*(t(i  ,j,k,nstp,itrc)-               &
!>   &                              t(i-1,j,k,nstp,itrc)))+             &
!>   &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
!>   &                    (0.75_r8*(tl_t(i  ,j,k,nrhs,itrc)-            &
!>   &                              tl_t(i-1,j,k,nrhs,itrc))+           &
!>   &                     0.25_r8*(tl_t(i  ,j,k,nstp,itrc)-            &
!>   &                              tl_t(i-1,j,k,nstp,itrc))))
!>
              adfac=cff*ad_FX(i,j)
              adfac1=adfac*(0.75_r8*(t(i  ,j,k,nrhs,itrc)-              &
     &                               t(i-1,j,k,nrhs,itrc))+             &
     &                      0.25_r8*(t(i  ,j,k,nstp,itrc)-              &
     &                               t(i-1,j,k,nstp,itrc)))
              adfac2=adfac*(Hz(i,j,k)+Hz(i-1,j,k))
              adfac3=adfac2*0.75_r8
              adfac4=adfac2*0.25_r8
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-adfac3
              ad_t(i  ,j,k,nrhs,itrc)=ad_t(i  ,j,k,nrhs,itrc)+adfac3
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)-adfac4
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+adfac4
              ad_FX(i,j)=0.0_r8
#elif defined TS_MIX_CLIMA
              IF (LtracerCLM(itrc,ng)) THEN
!>              tl_FX(i,j)=cff*                                         &
!>   &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
!>   &                      ((t(i  ,j,k,nrhs,itrc)-                     &
!>   &                        tclm(i  ,j,k,itrc))-                      &
!>   &                       (t(i-1,j,k,nrhs,itrc)-                     &
!>   &                        tclm(i-1,j,k,itrc)))+                     &
!>   &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
!>   &                      (tl_t(i  ,j,k,nrhs,itrc)-                   &
!>   &                       tl_t(i-1,j,k,nrhs,itrc)))
!>
                adfac=cff*ad_FX(i,j)
                adfac1=adfac*((t(i  ,j,k,nrhs,itrc)-                    &
     &                         tclm(i  ,j,k,itrc))-                     &
     &                        (t(i-1,j,k,nrhs,itrc)-                    &
     &                         tclm(i-1,j,k,itrc)))
                adfac2=adfac*(Hz(i,j,k)+Hz(i-1,j,k))
                ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
                ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
                ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-adfac2
                ad_t(i  ,j,k,nrhs,itrc)=ad_t(i  ,j,k,nrhs,itrc)+adfac2
                ad_FX(i,j)=0.0_r8
              ELSE
!>              tl_FX(i,j)=cff*                                         &
!>   &                     ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
!>   &                      (t(i  ,j,k,nrhs,itrc)-                      &
!>   &                       t(i-1,j,k,nrhs,itrc))+                     &
!>   &                      (Hz(i,j,k)+Hz(i-1,j,k))*                    &
!>   &                      (tl_t(i  ,j,k,nrhs,itrc)-                   &
!>   &                       tl_t(i-1,j,k,nrhs,itrc)))
!>
                adfac=cff*ad_FX(i,j)
                adfac1=adfac*(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
                adfac2=adfac*(Hz(i,j,k)+Hz(i-1,j,k))
                ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
                ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
                ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-adfac2
                ad_t(i  ,j,k,nrhs,itrc)=ad_t(i  ,j,k,nrhs,itrc)+adfac2
                ad_FX(i,j)=0.0_r8
              END IF
#else
!>            tl_FX(i,j)=cff*                                           &
!>   &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
!>   &                    (t(i  ,j,k,nrhs,itrc)-                        &
!>   &                     t(i-1,j,k,nrhs,itrc))+                       &
!>   &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
!>   &                    (tl_t(i  ,j,k,nrhs,itrc)-                     &
!>   &                     tl_t(i-1,j,k,nrhs,itrc)))
!>
              adfac=cff*ad_FX(i,j)
              adfac1=adfac*(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
              adfac2=adfac*(Hz(i,j,k)+Hz(i-1,j,k))
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-adfac2
              ad_t(i  ,j,k,nrhs,itrc)=ad_t(i  ,j,k,nrhs,itrc)+adfac2
              ad_FX(i,j)=0.0_r8
#endif
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE ad_t3dmix2_tile
