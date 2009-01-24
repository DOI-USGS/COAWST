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

      SUBROUTINE rp_t3dmix4 (ng, tile)
!
!svn $Id: rp_t3dmix4_s.h 694 2008-08-08 18:33:05Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes represneters tangent linear horizontal     !
!  biharmonic mixing of tracers along S-coordinate levels surfaces.    !
!                                                                      !
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iRPM, 27)
#endif
      CALL rp_t3dmix4_tile (ng, tile,                                   &
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
     &                      MIXING(ng) % diff4,                         &
#ifdef DIAGNOSTICS_TS
!!   &                      DIAGS(ng) % DiaTwrk,                        &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)
#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 27)
#endif
      RETURN
      END SUBROUTINE rp_t3dmix4
!
!***********************************************************************
      SUBROUTINE rp_t3dmix4_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nnew,                           &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
     &                            Hz, tl_Hz,                            &
     &                            pmon_u, pnom_v, pm, pn,               &
     &                            diff4,                                &
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
      real(r8), intent(in) :: diff4(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
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
      real(r8), intent(in) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
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
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
!  Compute horizontal tracer flux in the XI- and ETA-directions.
!
          DO j=J_RANGE
            DO i=I_RANGE+1
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*          &
     &            pmon_u(i,j)
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
              FX(i,j)=cff*(Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                (t(i  ,j,k,nrhs,itrc)-                            &
     &                 t(i-1,j,k,nrhs,itrc))
              tl_FX(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (t(i  ,j,k,nrhs,itrc)-                        &
     &                     t(i-1,j,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_t(i  ,j,k,nrhs,itrc)-                     &
     &                     tl_t(i-1,j,k,nrhs,itrc)))-                   &
#ifdef TL_IOMS
     &                   FX(i,j)
#endif
            END DO
          END DO
          DO j=J_RANGE+1
            DO i=I_RANGE
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*          &
                  pnom_v(i,j)
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
              FE(i,j)=cff*(Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                (t(i,j  ,k,nrhs,itrc)-                            &
     &                 t(i,j-1,k,nrhs,itrc))
              tl_FE(i,j)=cff*                                           &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (t(i,j  ,k,nrhs,itrc)-                        &
     &                     t(i,j-1,k,nrhs,itrc))+                       &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_t(i,j  ,k,nrhs,itrc)-                     &
     &                     tl_t(i,j-1,k,nrhs,itrc)))-                   &
#ifdef TL_IOMS
     &                   FE(i,j)
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
# ifndef EW_PERIODIC
          IF (WESTERN_EDGE) THEN
            DO j=J_RANGE
#  ifdef WESTERN_WALL
              LapT(Istr-1,j)=0.0_r8
              tl_LapT(Istr-1,j)=0.0_r8
#  else
              LapT(Istr-1,j)=LapT(Istr,j)
              tl_LapT(Istr-1,j)=tl_LapT(Istr,j)
#  endif
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=J_RANGE
#  ifdef EASTERN_WALL
              LapT(Iend+1,j)=0.0_r8
              tl_LapT(Iend+1,j)=0.0_r8
#  else
              LapT(Iend+1,j)=LapT(Iend,j)
              tl_LapT(Iend+1,j)=tl_LapT(Iend,j)
#  endif
            END DO
          END IF
# endif
# ifndef NS_PERIODIC
          IF (SOUTHERN_EDGE) THEN
            DO i=I_RANGE
#  ifdef SOUTHERN_WALL
              LapT(i,Jstr-1)=0.0_r8
              tl_LapT(i,Jstr-1)=0.0_r8
#  else
              LapT(i,Jstr-1)=LapT(i,Jstr)
              tl_LapT(i,Jstr-1)=tl_LapT(i,Jstr)
#  endif
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO i=I_RANGE
#  ifdef NORTHERN_WALL
              LapT(i,Jend+1)=0.0_r8
              tl_LapT(i,Jend+1)=0.0_r8
#  else
              LapT(i,Jend+1)=LapT(i,Jend)
              tl_LapT(i,Jend+1)=tl_LapT(i,Jend)
#  endif
            END DO
          END IF
# endif
!
!  Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*          &
     &            pmon_u(i,j)
!>            FX(i,j)=cff*                                              &
!>   &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
!>   &                (LapT(i,j)-LapT(i-1,j))
!>
              tl_FX(i,j)=cff*
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
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*          &
     &            pnom_v(i,j)
!>            FE(i,j)=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*      &
!>   &                (Hz(i,j,k)+Hz(i,j-1,k))*pnom_v(i,j)*              &
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
#ifdef TS_MPDATA_NOT_YET
              cff1=1.0_r8/Hz(i,j,k)
              tl_cff1=-cff1*cff1*tl_Hz(i,j,k)+                          &
# ifdef TL_IOMS
     &                2.0_r8*cff1
# endif
!>            t(i,j,k,3,itrc)=cff1*t(i,j,k,nnew,itrc)
!>
              tl_t(i,j,k,3,itrc)=tl_cff1*t(i,j,k,nnew,itrc)+            &
     &                           cff1*tl_t(i,j,k,nnew,itrc)-            &
# ifdef TL_IOMS
     &                           cff1*t(i,j,k,nnew,itrc)
# endif
#endif
#ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iThdif)=-cff
#endif
            END DO
          END DO
        END DO
      END DO
#undef I_RANGE
#undef J_RANGE
      RETURN
      END SUBROUTINE rp_t3dmix4_tile
