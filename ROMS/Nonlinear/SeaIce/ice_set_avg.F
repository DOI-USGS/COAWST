#include "cppdefs.h"
      MODULE ice_set_avg_mod
#if defined ICE_MODEL && defined SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module accumulates and computes output time-averaged ice       !
!  model fields. Due to synchronization, the time-averaged fields      !
!  are computed in delayed mode. All averages are accumulated at       !
!  the beggining of the next time-step.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC :: ice_set_avg
      PUBLIC :: ice_set_avg_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_set_avg (ng, tile, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 5, __LINE__, MyFile)
# endif
      CALL ice_set_avg_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       IOUT)

# ifdef PROFILE
      CALL wclock_off (ng, model, 5, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_set_avg
!
!***********************************************************************
      SUBROUTINE ice_set_avg_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Iout)
!***********************************************************************
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Iout
!
!  Local variable declarations.
!
      integer :: i, ifield, j, nf, ns
!
      real(r8) :: fac
!
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize ice model state variables.
!
        DO ns=1,nIceS
          IF (iSice(ns).gt.0) THEN
            ifield=iSice(ns)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_SAVG(ns,ng)%var(i,j)=ICE(ng)%Si(i,j,Iout,ns)
# ifdef MASKING
                  ICE_SAVG(ns,ng)%var(i,j)=ICE_SAVG(ns,ng)%var(i,j)*    &
     &                                     GRID(ng)%rmask_full(i,j)
# endif
                END DO
              END DO
            END IF
          END IF
        END DO
!
!  Initialize ice model internal variables.
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_FAVG(nf,ng)%var(i,j)=ICE(ng)%Fi(i,j,nf)
# ifdef MASKING
                  ICE_FAVG(nf,ng)%var(i,j)=ICE_FAVG(nf,ng)%var(i,j)*    &
     &                                     GRID(ng)%rmask_full(i,j)
# endif
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumulate ice model state variables.
!
        DO ns=1,nIceS
          IF (iSice(ns).gt.0) THEN
            ifield=iSice(ns)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_SAVG(ns,ng)%var(i,j)=ICE_SAVG(ns,ng)%var(i,j)+    &
# ifdef MASKING
     &                                     GRID(ng)%rmask_full(i,j)*    &
# endif
     &                                     ICE(ng)%Si(i,j,Iout,ns)
                END DO
              END DO
            END IF
          END IF
        END DO
!
!  Accumulate ice model internal variables.
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_FAVG(nf,ng)%var(i,j)=ICE_FAVG(nf,ng)%var(i,j)+    &
# ifdef MASKING
     &                                     GRID(ng)%rmask_full(i,j)*    &
# endif
     &                                     ICE(ng)%Fi(i,j,nf)
                END DO
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!  Notice that we need to apply periodic conditions, if any, since
!  the full I- and J-ranges are different.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
# ifdef WET_DRY
!  In wetting and drying, the sums are devided by the number of times
!  that each qrid point is wet.
# endif
!
# ifdef WET_DRY
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_avg(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_avg(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
# endif
!
!  Process ice model state variables.
!
        DO ns=1,nIceS
          IF (iSice(ns).gt.0) THEN
            ifield=iSice(ns)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_SAVG(ns,ng)%var(i,j)=rfac(i,j)*                   &
     &                                     ICE_SAVG(ns,ng)%var(i,j)
                END DO
              END DO
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  ICE_SAVG(ns,ng)%var)

# ifdef DISTRIBUTE
                CALL mp_exchange2d (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              ICE_SAVG(ns,ng)%var)
# endif
              END IF
            END IF
          END IF
        END DO
!
!  Process ice model internal variables.
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (Aout(ifield,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ICE_FAVG(nf,ng)%var(i,j)=rfac(i,j)*                   &
     &                                     ICE_FAVG(nf,ng)%var(i,j)
                END DO
              END DO
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  ICE_FAVG(ns,ng)%var)

# ifdef DISTRIBUTE
                CALL mp_exchange2d (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              ICE_FAVG(ns,ng)%var)
# endif
              END IF
            END IF
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE ice_set_avg_tile
#endif
      END MODULE ice_set_avg_mod
