#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE set_diags (ng, tile)
!
!svn $Id: set_diags.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  diagnostic fields.  Due to synchronization, the time-averaged       !
!  diagnostic fields are computed in delayed mode. All averages        !
!  are accumulated at the beginning of the next time-step.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_diags_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# ifdef SOLVE3D
     &                     kstp(ng),                                    &
# else
     &                     knew(ng),                                    &
# endif
     &                     nrhs(ng))
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_diags
!
!***********************************************************************
      SUBROUTINE set_diags_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kout, nrhs)
!***********************************************************************
!
      USE mod_param
# if defined DIAGNOSTICS_BIO
      USE mod_biology
# endif
      USE mod_diags
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
# endif
# ifdef WET_DRY
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
#  if defined WET_DRY || defined DIAGNOSTICS_BIO
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange4d
#  endif
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kout, nrhs
!
!  Local variable declarations.
!
      integer :: i, it, j, k
      integer :: idiag

      real(r8) :: fac

      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nDIA(ng).eq.0) RETURN

# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
!
!-----------------------------------------------------------------------
! Initialize time-averaged diagnostic arrays when appropriate.  Notice
! that fields are initilized twice during re-start.  However, the time-
! averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nDIA(ng)).eq.1)).or.                          &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN

#  ifdef WET_DRY
!
!  If wetting and drying, initialize time dependent counters for wet
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_dia(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_dia(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_dia(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_dia(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
#  endif
!
!  Initialize.
!
#  ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
#   ifdef WET_DRY
     &                      GRID(ng)%rmask_full(i,j)*                   &
#   endif
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef DIAGNOSTICS_UV
        DO idiag=1,NDM2d
          DO j=JstrR,JendR
            DO i=Istr,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=                              &
#   ifdef WET_DRY
     &                  GRID(ng)%umask_full(i,j)*                       &
#   endif
     &                  DIAGS(ng)%DiaU2wrk(i,j,idiag)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              DIAGS(ng)%DiaV2d(i,j,idiag)=                              &
#   ifdef WET_DRY
     &                  GRID(ng)%vmask_full(i,j)*                       &
#   endif
     &                  DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
        END DO
#   ifdef SOLVE3D
        DO idiag=1,NDM3d
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
#    ifdef WET_DRY
     &                    GRID(ng)%umask_full(i,j)*                     &
#    endif
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
              END DO
            END DO
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
#    ifdef WET_DRY
     &                    GRID(ng)%vmask_full(i,j)*                     &
#    endif
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
        END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsDIA(ng)) THEN

#  ifdef WET_DRY
!
!  If wetting and drying, accumulate wet points counters.
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_dia(i,j)=GRID(ng)%pmask_dia(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_dia(i,j)=GRID(ng)%rmask_dia(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_dia(i,j)=GRID(ng)%umask_dia(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_dia(i,j)=GRID(ng)%vmask_dia(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
#  endif
!
!  Accumulate diagnostic terms.
!
#  ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)+           &
#   ifdef WET_DRY
     &                      GRID(ng)%rmask_full(i,j)*                   &
#   endif
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef DIAGNOSTICS_UV
        DO idiag=1,NDM2d
          DO j=JstrR,JendR
            DO i=Istr,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)+  &
#   ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#   endif
     &                                    DIAGS(ng)%DiaU2wrk(i,j,idiag)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)+  &
#   ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#   endif
     &                                    DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
        END DO
#   ifdef SOLVE3D
        DO idiag=1,NDM3d
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaU3d(i,j,k,idiag)+                &
#    ifdef WET_DRY
     &                    GRID(ng)%umask_full(i,j)*                     &
#    endif
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
              END DO
            END DO
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaV3d(i,j,k,idiag)+                &
#    ifdef WET_DRY
     &                    GRID(ng)%vmask_full(i,j)*                     &
#    endif
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
        END DO
#   endif
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nDIA(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nDIA(ng).eq.1) THEN
            DIAtime(ng)=time(ng)
          ELSE
            DIAtime(ng)=DIAtime(ng)+REAL(nDIA(ng),r8)*dt(ng)
          END IF
        END IF
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
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_dia(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_dia(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_dia(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nDIA(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
# endif

# ifdef WET_DRY
!
!  If last time-step of average window, convert time dependent counters
!  for wet points to time-averaged Land/Sea masks (dry=0, wet=1) for
!  the current average window period. Notice that a grid point is wet
!  if the count is greater than zero for the current time average
!  window.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_dia(i,j)=MIN(1.0_r8, GRID(ng)%pmask_dia(i,j))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_dia(i,j)=MIN(1.0_r8, GRID(ng)%rmask_dia(i,j))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_dia(i,j)=MIN(1.0_r8, GRID(ng)%umask_dia(i,j))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_dia(i,j)=MIN(1.0_r8, GRID(ng)%vmask_dia(i,j))
          END DO
        END DO
!
!  Exchange boundary data.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            GRID(ng)%pmask_dia)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            GRID(ng)%rmask_dia)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            GRID(ng)%umask_dia)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            GRID(ng)%vmask_dia)
        END IF

#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      GRID(ng)%pmask_dia,                         &
     &                      GRID(ng)%rmask_dia,                         &
     &                      GRID(ng)%umask_dia,                         &
     &                      GRID(ng)%vmask_dia)
#  endif
# endif
!
!  Convert accumulated sums to time averages.
!
# ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=rfac(i,j)*           &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
# endif
# if defined DIAGNOSTICS_BIO
        DO idiag=1,NDbio2d
          IF (idiag.ne.ipCO2) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                DIAGS(ng)%DiaBio2d(i,j,idiag)=rfac(i,j)*                &
     &                    DIAGS(ng)%DiaBio2d(i,j,idiag)
              END DO
            END DO
          END IF
        END DO
        DO idiag=1,NDbio3d
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                DIAGS(ng)%DiaBio3d(i,j,k,idiag)=rfac(i,j)*              &
     &                    DIAGS(ng)%DiaBio3d(i,j,k,idiag)
              END DO
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
        DO idiag=1,NDM2d
          DO j=JstrR,JendR
            DO i=Istr,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=ufac(i,j)*                    &
     &                                    DIAGS(ng)%DiaU2d(i,j,idiag)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              DIAGS(ng)%DiaV2d(i,j,idiag)=vfac(i,j)*                    &
     &                                    DIAGS(ng)%DiaV2d(i,j,idiag)
            END DO
          END DO
        END DO
#  ifdef SOLVE3D
        DO idiag=1,NDM3d
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=ufac(i,j)*                &
     &                    DIAGS(ng)%DiaU3d(i,j,k,idiag)
              END DO
            END DO
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=vfac(i,j)*                &
     &                    DIAGS(ng)%DiaV3d(i,j,k,idiag)
              END DO
            END DO
          END DO
        END DO
#  endif
# endif

# if defined DIAGNOSTICS_TS  || defined DIAGNOSTICS_UV || \
     defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions for output purposes.
!-----------------------------------------------------------------------
!
#  ifdef DIAGNOSTICS_TS
!
!  3D tracer diagnostics.
!
        DO idiag=1,NDT
          DO it=1,NT(ng)
            CALL bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        DIAGS(ng)%DiaTrc(:,:,:,it,idiag))
          END DO
#   ifdef DISTRIBUTE
          CALL mp_exchange4d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),  &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        DIAGS(ng)%DiaTrc(:,:,:,:,idiag))
#   endif
        END DO
#  endif
#  ifdef DIAGNOSTICS_UV
!
!  2D momentum diagnostics.
!
        DO idiag=1,NDM2d
          CALL bc_u2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      DIAGS(ng)%DiaU2d(:,:,idiag))
          CALL bc_v2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      DIAGS(ng)%DiaV2d(:,:,idiag))
        END DO
#   ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, NDM2d,               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      DIAGS(ng)%DiaU2d,                           &
     &                      DIAGS(ng)%DiaV2d)
#   endif
#   ifdef SOLVE3D
!
!  3D momentum diagnostics.
!
        DO idiag=1,NDM3d
          CALL bc_u3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      DIAGS(ng)%DiaU3d(:,:,:,idiag))
          CALL bc_v3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      DIAGS(ng)%DiaV3d(:,:,:,idiag))
        END DO
#    ifdef DISTRIBUTE
        CALL mp_exchange4d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NDM3d,     &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      DIAGS(ng)%DiaU3d,                           &
     &                      DIAGS(ng)%DiaV3d)
#    endif
#   endif
#  endif
#  if defined DIAGNOSTICS_BIO
!
!  Biological terms 2D diagnostics.
!
        DO idiag=1,NDbio2d
          CALL bc_r2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      DIAGS(ng)%DiaBio2d(:,:,idiag))
        END DO
#   ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, NDbio2d,             &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      DIAGS(ng)%DiaBio2d)
#   endif
!
!  Biological terms 3D diagnostics.
!
        DO idiag=1,NDbio3d
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      DIAGS(ng)%DiaBio3d(:,:,:,idiag))
        END DO
#   ifdef DISTRIBUTE
        CALL mp_exchange4d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NDbio3d,   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      DIAGS(ng)%DiaBio3d)
#   endif
#  endif
# endif
      END IF
      RETURN
      END SUBROUTINE set_diags_tile
#else
      SUBROUTINE set_diags
      RETURN
      END SUBROUTINE set_diags
#endif
