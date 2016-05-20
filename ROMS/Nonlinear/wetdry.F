#include "cppdefs.h"
      MODULE wetdry_mod
#ifdef WET_DRY
!
!svn $Id: exchange_2d.F 503 2008-01-10 00:11:51Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes the wet/dry masking arrays.                   !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wetdry (ng, tile, Tindex, Linitialize)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: Linitialize
      integer, intent(in) :: ng, tile, Tindex
!
!  Local variable declarations.
!
# include "tile.h"
!
      IF (Linitialize) THEN
        CALL wetdry_ini_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        GRID(ng)%pmask,                           &
     &                        GRID(ng)%rmask,                           &
     &                        GRID(ng)%umask,                           &
     &                        GRID(ng)%vmask,                           &
# endif
     &                        GRID(ng)%h,                               &
     &                        OCEAN(ng)%zeta(:,:,Tindex),               &
# ifdef SOLVE3D
     &                        OCEAN(ng)%ubar(:,:,Tindex),               &
     &                        OCEAN(ng)%vbar(:,:,Tindex),               &
     &                        GRID(ng)%umask_diff,                      &
     &                        GRID(ng)%vmask_diff,                      &
# endif
     &                        GRID(ng)%pmask_wet,                       &
     &                        GRID(ng)%pmask_full,                      &
     &                        GRID(ng)%rmask_wet,                       &
     &                        GRID(ng)%rmask_full,                      &
     &                        GRID(ng)%umask_wet,                       &
     &                        GRID(ng)%umask_full,                      &
     &                        GRID(ng)%vmask_wet,                       &
     &                        GRID(ng)%vmask_full)
      ELSE
        CALL wetdry_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    GRID(ng)%pmask,                               &
     &                    GRID(ng)%rmask,                               &
     &                    GRID(ng)%umask,                               &
     &                    GRID(ng)%vmask,                               &
# endif
     &                    GRID(ng)%h,                                   &
     &                    OCEAN(ng)%zeta(:,:,Tindex),                   &
# ifdef SOLVE3D
     &                    COUPLING(ng)%DU_avg1,                         &
     &                    COUPLING(ng)%DV_avg1,                         &
     &                    GRID(ng)%rmask_wet_avg,                       &
     &                    GRID(ng)%umask_diff,                          &
     &                    GRID(ng)%vmask_diff,                          &
# endif
     &                    GRID(ng)%pmask_wet,                           &
     &                    GRID(ng)%pmask_full,                          &
     &                    GRID(ng)%rmask_wet,                           &
     &                    GRID(ng)%rmask_full,                          &
     &                    GRID(ng)%umask_wet,                           &
     &                    GRID(ng)%umask_full,                          &
     &                    GRID(ng)%vmask_wet,                           &
     &                    GRID(ng)%vmask_full)
      END IF

      RETURN
      END SUBROUTINE wetdry
!
!***********************************************************************
      SUBROUTINE wetdry_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        pmask, rmask, umask, vmask,               &
# endif
     &                        h, zeta,                                  &
# ifdef SOLVE3D
     &                        DU_avg1, DV_avg1,                         &
     &                        rmask_wet_avg,                            &
     &                        umask_diff, vmask_diff,                   &
# endif
     &                        pmask_wet, pmask_full,                    &
     &                        rmask_wet, rmask_full,                    &
     &                        umask_wet, umask_full,                    &
     &                        vmask_wet, vmask_full)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:,LBj:)
      real(r8), intent(out) :: umask_diff(LBi:,LBj:)
      real(r8), intent(out) :: vmask_diff(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: pmask_full(LBi:,LBj:)
      real(r8), intent(out) :: rmask_full(LBi:,LBj:)
      real(r8), intent(out) :: umask_full(LBi:,LBj:)
      real(r8), intent(out) :: vmask_full(LBi:,LBj:)

      real(r8), intent(out) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)

# else

      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_diff(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_diff(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: pmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_full(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j

      real(r8) :: cff
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
! Compute local mask at RHO-points in terms of free-surface. The local
! array allows computations in the same parallel region.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          wetdry(i,j)=1.0_r8
# ifdef MASKING
          wetdry(i,j)=wetdry(i,j)*rmask(i,j)
# endif
          IF ((zeta(i,j)+h(i,j)).le.(Dcrit(ng)+eps)) THEN
            wetdry(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
! Compute wet/dry mask arrays.
!
      IF (iif(ng).le.nfast(ng)) THEN
        CALL wetdry_mask_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         wetdry,                                  &
     &                         pmask_wet, rmask_wet,                    &
     &                         umask_wet, vmask_wet)
      END IF

# ifdef SOLVE3D
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
      IF (iif(ng).le.nfast(ng)) THEN
        IF (PREDICTOR_2D_STEP(ng).and.(FIRST_2D_STEP)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=wetdry(i,j)
            END DO
          END DO
        ELSE
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=rmask_wet_avg(i,j)+wetdry(i,j)
            END DO
          END DO
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_wet_avg)
        END IF

#  ifdef DISTRIBUTE
!
        CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      rmask_wet_avg)
#  endif
!
!  If done witn 2D fast time-stepping, scale time-averaged mask by
!  the inverse of fast time-step (2*nfast).
!
      ELSE
        cff=1.0_r8/REAL(2*nfast(ng),r8)
        DO j=Jstr-1,JendR                      ! executed in a different
          DO i=Istr-1,IendR                    ! paralllel region
            wetdry(i,j)=AINT(rmask_wet_avg(i,j)*cff)
          END DO
        END DO
!
! Compute time averaged wet/dry mask arrays.
!
        CALL wetdry_avg_mask_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             DU_avg1, DV_avg1,                    &
     &                             wetdry,                              &
     &                             pmask_wet, rmask_wet,                &
     &                             umask_wet, vmask_wet,                &
     &                             umask_diff, vmask_diff)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set masks full time-dependent masks.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      IF (iif(ng).gt.nfast(ng)) THEN
# else
      IF (iif(ng).eq.nfast(ng)) THEN
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_full(i,j)=rmask_wet(i,j)*rmask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_full(i,j)=MAX(pmask_wet(i,j)*pmask(i,j), 2.0_r8)
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_full(i,j)=umask_wet(i,j)*umask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_full(i,j)=vmask_wet(i,j)*vmask(i,j)
          END DO
        END DO
!
!  Insure that masks at mass point source locations are set to water
!  to avoid writting output with FillValue at those locations.
!
        IF (LuvSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                umask_full(i,j)=1.0_r8
              ELSE
                vmask_full(i,j)=1.0_r8
              END IF
            END IF
          END DO
        END IF
!
!  Exchange boundary data.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pmask_full)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_full)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            umask_full)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vmask_full)
        END IF

# ifdef DISTRIBUTE
!
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      pmask_full, rmask_full,                     &
     &                      umask_full, vmask_full)
# endif
      END IF

      RETURN
      END SUBROUTINE wetdry_tile
!
!***********************************************************************
      SUBROUTINE wetdry_ini_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef MASKING
     &                            pmask, rmask, umask, vmask,           &
# endif
     &                            h, zeta,                              &
# ifdef SOLVE3D
     &                            ubar, vbar,                           &
     &                            umask_diff, vmask_diff,               &
# endif
     &                            pmask_wet, pmask_full,                &
     &                            rmask_wet, rmask_full,                &
     &                            umask_wet, umask_full,                &
     &                            vmask_wet, vmask_full)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ubar(LBi:,LBj:)
      real(r8), intent(in) :: vbar(LBi:,LBj:)
      real(r8), intent(out) :: umask_diff(LBi:,LBj:)
      real(r8), intent(out) :: vmask_diff(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: pmask_full(LBi:,LBj:)
      real(r8), intent(out) :: rmask_full(LBi:,LBj:)
      real(r8), intent(out) :: umask_full(LBi:,LBj:)
      real(r8), intent(out) :: vmask_full(LBi:,LBj:)

      real(r8), intent(out) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)

# else

      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_diff(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_diff(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: pmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_full(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j

      real(r8) :: cff
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
! Compute local mask at RHO-points in terms of free-surface. The local
! array allows computations in the same parallel region.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          wetdry(i,j)=1.0_r8
# ifdef MASKING
          wetdry(i,j)=wetdry(i,j)*rmask(i,j)
# endif
          IF ((zeta(i,j)+h(i,j)).le.(Dcrit(ng)+eps)) THEN
            wetdry(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
! Compute initial wet/dry mask arrays.
!
# ifdef SOLVE3D
      CALL wetdry_avg_mask_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ubar, vbar,                            &
     &                           wetdry,                                &
     &                           pmask_wet, rmask_wet,                  &
     &                           umask_wet, vmask_wet,                  &
     &                           umask_diff, vmask_diff)
# else
      CALL wetdry_mask_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       wetdry,                                    &
     &                       pmask_wet, rmask_wet,                      &
     &                       umask_wet, vmask_wet)
# endif
!
!-----------------------------------------------------------------------
!  Set masks full time-dependent masks.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_full(i,j)=rmask_wet(i,j)*rmask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_full(i,j)=MAX(pmask_wet(i,j)*pmask(i,j), 2.0_r8)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_full(i,j)=umask_wet(i,j)*umask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_full(i,j)=vmask_wet(i,j)*vmask(i,j)
        END DO
      END DO
!
!  Insure that masks at mass point source locations are set to water
!  to avoid writting output with FillValue at those locations.
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              umask_full(i,j)=1.0_r8
            ELSE
              vmask_full(i,j)=1.0_r8
            END IF
          END IF
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_full)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_full)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_full)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_full)
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_full, rmask_full,                       &
     &                    umask_full, vmask_full)
# endif

      RETURN
      END SUBROUTINE wetdry_ini_tile
!
!***********************************************************************
      SUBROUTINE wetdry_mask_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             wetdry,                              &
     &                             pmask_wet, rmask_wet,                &
     &                             umask_wet, vmask_wet)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: wetdry(IminS:,JminS:)
      real(r8), intent(out) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)
# else
      real(r8), intent(in)  :: wetdry(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute wet/dry masks in terms of RHO-points mask (wetdry).
!-----------------------------------------------------------------------
!
!  Wet/dry mask at RHO-points.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_wet(i,j)=wetdry(i,j)
        END DO
      END DO
!
!  Wet/dry mask at U-points.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_wet(i,j)=wetdry(i-1,j)+wetdry(i,j)
          IF (umask_wet(i,j).eq.1.0_r8) THEN
            umask_wet(i,j)=wetdry(i-1,j)-wetdry(i,j)
          END IF
        END DO
      END DO
!
!  Wet/dry mask at V-points.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_wet(i,j)=wetdry(i,j-1)+wetdry(i,j)
          IF (vmask_wet(i,j).eq.1.0_r8) THEN
            vmask_wet(i,j)=wetdry(i,j-1)-wetdry(i,j)
          END IF
        END DO
      END DO
!
!  Wet/dry mask at PSI-points.
!
      cff1=1.0_r8
      cff2=2.0_r8
      DO j=Jstr,JendR
        DO i=Istr,IendR
          IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                          &
     &        (wetdry(i  ,j  ).gt.0.5_r8).and.                          &
     &        (wetdry(i-1,j-1).gt.0.5_r8).and.                          &
     &        (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=1.0_r8
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE
            pmask_wet(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_wet)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_wet)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_wet)
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_wet, rmask_wet, umask_wet, vmask_wet)
# endif

      RETURN
      END SUBROUTINE wetdry_mask_tile

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE wetdry_avg_mask_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 DU_avg1, DV_avg1,                &
     &                                 wetdry,                          &
     &                                 pmask_wet, rmask_wet,            &
     &                                 umask_wet, vmask_wet,            &
     &                                 umask_diff, vmask_diff)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: wetdry(IminS:,JminS:)
      real(r8), intent(in)  :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in)  :: DV_avg1(LBi:,LBj:)
      real(r8), intent(out) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_diff(LBi:,LBj:)
      real(r8), intent(out) :: vmask_diff(LBi:,LBj:)
#  else
      real(r8), intent(in)  :: wetdry(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in)  :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_diff(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_diff(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2, cff5, cff6

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time-averaged wet/dry masks in terms of RHO-points mask
!  (wetdry).
!-----------------------------------------------------------------------
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_wet(i,j)=wetdry(i,j)
        END DO
      END DO
!
!  Wet/dry mask at U-points, averaged over all fast time-steps.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff1=wetdry(i-1,j)+wetdry(i,j)
          umask_diff(i,j) = 1.0_r8
          IF (cff1.lt.2.0_r8) THEN
            umask_diff(i,j) = 0.0_r8
          END IF
          IF (cff1.eq.1.0_r8) THEN
            cff1=wetdry(i-1,j)-wetdry(i,j)
          END IF
          cff5=ABS(ABS(cff1)-1.0_r8)
          cff6=0.5_r8+DSIGN(0.5_r8,DU_avg1(i,j))*cff1
          umask_wet(i,j)=0.5_r8*cff1*cff5+cff6*(1.0_r8-cff5)
!                                                       catch lone ponds
          IF (DU_avg1(i,j).eq.0.0_r8) THEN
            IF ((wetdry(i-1,j)+wetdry(i,j)).le.1.0_r8) THEN
              umask_wet(i,j)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Wet/dry mask at V-points, averaged over all fast time-steps.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff1=wetdry(i,j-1)+wetdry(i,j)
          vmask_diff(i,j) = 1.0_r8
          IF (cff1.lt.2.0_r8) THEN
            vmask_diff(i,j) = 0.0_r8
          END IF
          IF (cff1.eq.1.0_r8) THEN
            cff1=wetdry(i,j-1)-wetdry(i,j)
          END IF
          cff5=ABS(ABS(cff1)-1.0_r8)
          cff6=0.5_r8+DSIGN(0.5_r8,DV_avg1(i,j))*cff1
          vmask_wet(i,j)=0.5_r8*cff1*cff5+cff6*(1.0_r8-cff5)
!                                                       catch lone ponds
          IF (DV_avg1(i,j).eq.0.0_r8) THEN
            IF ((wetdry(i,j-1)+wetdry(i,j)).le.1.0_r8) THEN
              vmask_wet(i,j)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Wet/dry mask at PSI-points.
!
      cff1=1.0_r8
      cff2=2.0_r8
      DO j=Jstr,JendR
        DO i=Istr,IendR
          IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                          &
     &        (wetdry(i  ,j  ).gt.0.5_r8).and.                          &
     &        (wetdry(i-1,j-1).gt.0.5_r8).and.                          &
     &        (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=1.0_r8
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff1
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).gt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).lt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE IF ((wetdry(i-1,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i  ,j  ).lt.0.5_r8).and.                     &
     &             (wetdry(i-1,j-1).gt.0.5_r8).and.                     &
     &             (wetdry(i  ,j-1).gt.0.5_r8)) THEN
            pmask_wet(i,j)=cff2
          ELSE
            pmask_wet(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_wet)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_wet)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_diff)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_diff)
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_wet, rmask_wet, umask_wet, vmask_wet)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    umask_diff, vmask_diff)
#  endif

      RETURN
      END SUBROUTINE wetdry_avg_mask_tile
# endif
#endif
      END MODULE wetdry_mod
