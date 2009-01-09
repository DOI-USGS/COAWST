#include "cppdefs.h"
#if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
      SUBROUTINE set_diags (ng, tile)
!
!svn $Id: set_diags.F 1598 2008-07-25 03:29:09Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
     &                     LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_diags
!
!***********************************************************************
      SUBROUTINE set_diags_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_diags
      USE mod_grid
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, it, j, k
      integer :: idiag

      real(r8) :: fac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nDIA(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
! Initialize time-averaged diagnostic arrays when appropriate.  Notice
! that fields are initilized twice during re-start.  However, the time-
! averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nDIA(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
# ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#  ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
#  endif
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsDIA(ng)) THEN
# ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)+           &
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#  ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaU3d(i,j,k,idiag)+                &
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaV3d(i,j,k,idiag)+                &
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions and land-mask
!  for output purposes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsDIA(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nDIA(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        IF (SOUTH_WEST_TEST) THEN
          DIAtime(ng)=DIAtime(ng)+REAL(nDIA(ng),r8)*dt(ng)
        END IF
# ifdef DIAGNOSTICS_TS
!
!  Apply periodic or gradient boundary conditions for output purposes.
!
        DO idiag=1,NDT
          DO it=1,NT(ng)
            CALL bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        DIAGS(ng)%DiaTrc(:,:,:,it,idiag))
#  ifdef MASKING
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)*           &
     &                      GRID(ng)%rmask(i,j)
                END DO
              END DO
            END DO
#  endif
          END DO
#  ifdef DISTRIBUTE
          CALL mp_exchange4d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),  &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        DIAGS(ng)%DiaTrc(:,:,:,:,idiag))
#  endif
        END DO
# endif
# ifdef DIAGNOSTICS_UV
!
!  Apply periodic or gradient boundary conditions for output purposes.
!
        DO idiag=1,NDM2d
          CALL bc_u2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      DIAGS(ng)%DiaU2d(:,:,idiag))
          CALL bc_v2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      DIAGS(ng)%DiaV2d(:,:,idiag))
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, NDM2d,               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%DiaU2d,                           &
     &                      DIAGS(ng)%DiaV2d)
#  endif
#  ifdef SOLVE3D
        DO idiag=1,NDM3d
          CALL bc_u3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      DIAGS(ng)%DiaU3d(:,:,:,idiag))
          CALL bc_v3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      DIAGS(ng)%DiaV3d(:,:,:,idiag))
        END DO
#   ifdef DISTRIBUTE
        CALL mp_exchange4d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NDM3d,     &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%DiaU3d,                           &
     &                      DIAGS(ng)%DiaV3d)
#   endif
#  endif
#  ifdef MASKING
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)*  &
     &                                    GRID(ng)%umask(i,j)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)*  &
     &                                    GRID(ng)%vmask(i,j)
            END DO
          END DO
#   ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                                 DIAGS(ng)%DiaU3d(i,j,k,idiag)*   &
     &                                 GRID(ng)%umask(i,j)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                                 DIAGS(ng)%DiaV3d(i,j,k,idiag)*   &
     &                                 GRID(ng)%vmask(i,j)
              END DO
            END DO
          END DO
#   endif
        END DO
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
