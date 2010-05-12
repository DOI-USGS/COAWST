#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE set_diags (ng, tile)
!
!svn $Id: set_diags.F 1598 2008-07-25 03:29:09Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      USE mod_diags
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
#  ifdef SOLVE3D
      USE exchange_3d_mod, ONLY : exchange_p3d_tile
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange4d
#  endif
# endif
      USE vorticity_mod, ONLY : vorticity_tile
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
# ifdef SOLVE3D
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nDIA(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity diagnostics.
!-----------------------------------------------------------------------
!
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef SOLVE3D
     &                       kout, nrhs,                                &
# else
     &                       kout,                                      &
# endif
# ifdef MASKING
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % fomn,                           &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef SOLVE3D
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % pden,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
# endif
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
# ifdef SOLVE3D
     &                       potvor, relvor,                            &
# endif
                             potvor_bar, relvor_bar)

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
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO j=Jstr,JendR
          DO i=Istr,IendR
            DIAGS(ng)%pvor_bar(i,j)=potvor_bar(i,j)
            DIAGS(ng)%rvor_bar(i,j)=relvor_bar(i,j)
          END DO
        END DO
#  ifdef SOLVE3D
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=Istr,IendR
              DIAGS(ng)%pvor(i,j,k)=potvor(i,j,k)
              DIAGS(ng)%rvor(i,j,k)=relvor(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef DIAGNOSTICS_TS
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
#  endif
#  ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#   ifdef SOLVE3D
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
#   endif
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsDIA(ng)) THEN
        DO j=Jstr,JendR
          DO i=Istr,IendR
            DIAGS(ng)%pvor_bar(i,j)=DIAGS(ng)%pvor_bar(i,j)+            &
     &                              potvor_bar(i,j)
            DIAGS(ng)%rvor_bar(i,j)=DIAGS(ng)%rvor_bar(i,j)+            &
     &                              relvor_bar(i,j)
          END DO
        END DO
#  ifdef SOLVE3D
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=Istr,IendR
              DIAGS(ng)%pvor(i,j,k)=DIAGS(ng)%pvor(i,j,k)+              &
     &                              potvor(i,j,k)
              DIAGS(ng)%rvor(i,j,k)=DIAGS(ng)%rvor(i,j,k)+              &
     &                              relvor(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef DIAGNOSTICS_TS
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
#  endif
#  ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#   ifdef SOLVE3D
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
#   endif
        END DO
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set diagnotics time.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsDIA(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nDIA(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        IF (SOUTH_WEST_TEST) THEN
          DIAtime(ng)=DIAtime(ng)+REAL(nDIA(ng),r8)*dt(ng)
        END IF

# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions and land-mask
!  for output purposes.
!-----------------------------------------------------------------------
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          DIAGS(ng)%pvor_bar)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          DIAGS(ng)%rvor_bar)
#  endif
#  ifdef MASKING
        DO j=Jstr,JendR
          DO i=Istr,IendR
            DIAGS(ng)%pvor_bar(i,j)=DIAGS(ng)%pvor_bar(i,j)*            &
     &                              GRID(ng)%pmask(i,j)
            DIAGS(ng)%rvor_bar(i,j)=DIAGS(ng)%rvor_bar(i,j)*            &
     &                              GRID(ng)%pmask(i,j)
          END DO
        END DO
#  endif
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%pvor_bar,                         &
     &                      DIAGS(ng)%rvor_bar)
#  endif
#  ifdef SOLVE3D
#   if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_p3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          DIAGS(ng)%pvor)
        CALL exchange_p3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          DIAGS(ng)%rvor)
#   endif
#   ifdef MASKING
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=Istr,IendR
              DIAGS(ng)%pvor(i,j,k)=DIAGS(ng)%pvor(i,j,k)*              &
     &                              GRID(ng)%pmask(i,j)
              DIAGS(ng)%rvor(i,j,k)=DIAGS(ng)%rvor(i,j,k)*              &
     &                              GRID(ng)%pmask(i,j)
            END DO
          END DO
        END DO
#   endif
#   ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%pvor,                             &
     &                      DIAGS(ng)%rvor)
#   endif
#  endif
#  ifdef DIAGNOSTICS_TS
!
!  3D diagnostics.
!
        DO idiag=1,NDT
          DO it=1,NT(ng)
            CALL bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        DIAGS(ng)%DiaTrc(:,:,:,it,idiag))
#   ifdef MASKING
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)*           &
     &                      GRID(ng)%rmask(i,j)
                END DO
              END DO
            END DO
#   endif
          END DO
#   ifdef DISTRIBUTE
          CALL mp_exchange4d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),  &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        DIAGS(ng)%DiaTrc(:,:,:,:,idiag))
#   endif
        END DO
#  endif
#  ifdef DIAGNOSTICS_UV
!
!  2D diagnostics.
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
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%DiaU2d,                           &
     &                      DIAGS(ng)%DiaV2d)
#   endif
#   ifdef SOLVE3D
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
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      DIAGS(ng)%DiaU3d,                           &
     &                      DIAGS(ng)%DiaV3d)
#    endif
#   endif
#   ifdef MASKING
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)*  &
     &                                    GRID(ng)%umask(i,j)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)*  &
     &                                    GRID(ng)%vmask(i,j)
            END DO
          END DO
#    ifdef SOLVE3D
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
#    endif
        END DO
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
