#include "cppdefs.h"
#ifdef FILTERED
      SUBROUTINE set_filter(ng, tile)
!
!=======================================================================
!  Copyright (c) 2002-2014 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beginning of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
# ifdef WET_DRY
      USE mod_grid
# endif
# ifdef WET_DRY
      USE set_masks_mod, ONLY : set_avg_masks
# endif
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5, __LINE__, __FILE__)
# endif
      CALL set_filter_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
# ifdef ICE_MODEL
     &                      IOUT, IUOUT, IEOUT,                         &
# endif
     &                      KOUT, NOUT)

# ifdef WET_DRY
      CALL set_avg_masks (ng, tile, iNLM,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % pmask_avg,                         &
     &                    GRID(ng) % rmask_avg,                         &
     &                    GRID(ng) % umask_avg,                         &
     &                    GRID(ng) % vmask_avg)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE set_filter
!
!***********************************************************************
      SUBROUTINE set_filter_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef ICE_MODEL
     &                            Iout, Iuout, Ieout,                   &
# endif
     &                            Kout, Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_filter
      USE mod_forces
      USE mod_grid
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
      USE mod_parallel
      USE uv_rotate_mod, ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod, ONLY : uv_rotate3d
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
# ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
# endif
      integer, intent(in) :: Kout, Nout
!
!  Local variable declarations.
!
      integer  :: i, itrc, j, k, iw, f
      real(r8) :: fac
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
# ifdef SOLVE3D
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
! Initialize time-filtered arrays when appropriate.  They aren't all
! initialized at beginning of run - skip right to summing, which is
! OK because it was initialized in init_arrays
!-----------------------------------------------------------------------
!
!
!  handle tracers other than T,S without filtering
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    (iic(ng).eq.ntstart(ng))) THEN
# ifdef WET_DRY
!
!  If wetting and drying, initialize time dependent counters for wet
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
        DO itrc=NAT+1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                         &
     &                    OCEAN(ng)%t(i,j,k,Nout,itrc)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                         &
     &                              AVERAGE(ng)%avgt(i,j,k,itrc)*       &
     &                              GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO
# ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=OCEAN(ng)%PONsed(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgPONsed(i,j)=AVERAGE(ng)%avgPONsed(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=OCEAN(ng)%OPALsed(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALsed(i,j)*  &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=OCEAN(ng)%DENITsed(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDENITsed(i,j)=AVERAGE(ng)%avgDENITsed(i,j)*&
     &                                     GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=OCEAN(ng)%PON_burial(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgPONbur(i,j)=AVERAGE(ng)%avgPONbur(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALbur(i,j)=OCEAN(ng)%OPAL_burial(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgOPALbur(i,j)=AVERAGE(ng)%avgOPALbur(i,j)*  &
     &                                    GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF
# endif
# ifdef PRIMARY_PROD
        IF (Aout(idNPP,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgNPP(i,j)=OCEAN(ng)%Bio_NPP(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgNPP(i,j)=AVERAGE(ng)%avgNPP(i,j)*          &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif

!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
# ifdef WET_DRY
!
!  If wetting and drying, accumulate wet points counters.
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg(i,j)=GRID(ng)%pmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg(i,j)=GRID(ng)%rmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg(i,j)=GRID(ng)%umask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg(i,j)=GRID(ng)%vmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
        DO itrc=NAT+1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                         &
     &                      AVERAGE(ng)%avgt(i,j,k,itrc)+               &
#  ifdef WET_DRY
     &                                       GRID(ng)%rmask_full(i,j)*  &
#  endif
     &                      OCEAN(ng)%t(i,j,k,Nout,itrc)
                END DO
              END DO
            END DO
          END IF
        END DO
#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=AVERAGE(ng)%avgPONsed(i,j)+    &
# ifdef WET_DRY
     &                                   GRID(ng)%rmask_full(i,j)*      &
# endif
     &                               OCEAN(ng)%PONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALsed(i,j)+  &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                               OCEAN(ng)%OPALsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=                             &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &           AVERAGE(ng)%avgDENITsed(i,j)+OCEAN(ng)%DENITsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=AVERAGE(ng)%avgPONbur(i,j)+    &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                               OCEAN(ng)%PON_burial(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALbur(i,j)+  &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                               OCEAN(ng)%OPAL_burial(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (Aout(idNPP,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgNPP(i,j)=AVERAGE(ng)%avgNPP(i,j)+          &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 OCEAN(ng)%Bio_NPP(i,j)
            END DO
          END DO
        END IF
#  endif
      END IF
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)))) THEN
# ifdef WET_DRY
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=fac
          END DO
        END DO
# endif
        DO itrc=NAT+1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                         &
     &                      rfac(i,j)*AVERAGE(ng)%avgt(i,j,k,itrc)
                END DO
              END DO
            END DO
          END IF
        END DO
#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=rfac(i,j)*                     &
     &                          AVERAGE(ng)%avgPONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=rfac(i,j)*                    &
     &                          AVERAGE(ng)%avgOPALsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=rfac(i,j)*                   &
     &                          AVERAGE(ng)%avgDENITsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=rfac(i,j)*                     &
     &                        AVERAGE(ng)%avgPONbur(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALbur(i,j)=rfac(i,j)*                    &
     &                        AVERAGE(ng)%avgOPALbur(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (Aout(idNPP,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgNPP(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgNPP(i,j)
            END DO
          END DO
        END IF
#  endif
      END IF
!     IF (Master) THEN
!       print *, 'SET_FILTER  nw:',nw,'dtfil',dtfil,'dt',INT(dt(ng))
!       print *, 'SET_FILTER  nfil:',nfil,'navg',navg,'nfil/nw',nfil/nw
!       print *, 'SET_FILTER  iic-1:',iic-1
!       print *, 'SET_FILTER  ntstart:',ntstart,'nrrec',nrrec
!     END IF
      IF (Aout(id2dPV,ng).or.Aout(id2dRV,ng).or.                        &
     &    Aout(id3dPV,ng).or.Aout(id3dRV,ng)) THEN
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef SOLVE3D
     &                       Kout, Nout,                                &
# else
     &                       Kout,                                      &
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
     &                       potvor_bar, relvor_bar)
        END IF
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &     (MOD(iic(ng)-ntstart(ng),nfil/nw).eq.0).and.                 &
     &    ((iic(ng).ne.ntstart(ng)))) THEN
!    &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN

        DO f = 1,FILN
          IF (fcount(f).eq.0) THEN
            iw=1
            IF (Master) THEN
              print *, 'SET_FILTER initialize sum', f, '    iw:', iw
            END IF
            IF (Aout(idFsur,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filzeta(i,j,f)=                            &
     &                                   OCEAN(ng)%zeta(i,j,Kout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,Iend+1
                  FILTER(ng)%filu2d (i,j,f)=                            &
     &                                   OCEAN(ng)%ubar(i,j,Kout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
              DO j=Jstr,Jend+1
                DO i=IstrR,IendR
                  FILTER(ng)%filv2d (i,j,f)=                            &
     &                                   OCEAN(ng)%vbar(i,j,Kout)*fw(iw)
                END DO
              END DO
            END IF
# ifdef SOLVE3D
            IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
              DO k=1,N(ng)
                DO j=JstrR,JendR
                  DO i=Istr,Iend+1
                    FILTER(ng)%filu3d(i,j,k,f)=                         &
     &                                   OCEAN(ng)%u(i,j,k,Nout)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
                  DO i=IstrR,IendR
                    FILTER(ng)%filv3d(i,j,k,f)=                         &
     &                                   OCEAN(ng)%v(i,j,k,Nout)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idDano,ng)) THEN
              DO k=1,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filrho(i,j,k,f)=                         &
     &                                   OCEAN(ng)%rho(i,j,k)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            DO itrc=1,NAT
              IF (Aout(idTvar(itrc),ng)) THEN
                DO k=1,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filt(i,j,k,itrc,f)=                    &
     &                         OCEAN(ng)%t(i,j,k,Nout,itrc)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
            END DO
            IF (Aout(idOvel,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filomega(i,j,k,f)=                       &
     &                                 OCEAN(ng)%W(i,j,k)*fw(iw)*       &
     &                                 GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idWvel,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filw3d(i,j,k,f)=                         &
     &                                 OCEAN(ng)%wvel(i,j,k)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
#  ifdef LMD_SKPP
            IF (Aout(idHsbl,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhsbl(i,j,f)=MIXING(ng)%hsbl(i,j)*fw(iw)
                END DO
              END DO
            END IF
#  endif
#  ifdef LMD_BKPP
            IF (Aout(idHbbl,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhbbl(i,j,f)=MIXING(ng)%hbbl(i,j)*fw(iw)
                END DO
              END DO
            END IF
#  endif
            IF (Aout(idVvis,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKv(i,j,k,f)=                         &
     &                                      MIXING(ng)%Akv(i,j,k)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idTdif,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKt(i,j,k,f)=                         &
     &                                MIXING(ng)%Akt(i,j,k,itemp)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idSdif,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKs(i,j,k,f)=                         &
     &                                MIXING(ng)%Akt(i,j,k,isalt)*fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idTsur(itemp),ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filstf(i,j,f)=                             &
     &                                FORCES(ng)%stflx(i,j,itemp)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idTsur(isalt),ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filswf(i,j,f)=                             &
     &                                FORCES(ng)%stflx(i,j,isalt)*fw(iw)
                END DO
              END DO
            END IF
#  ifdef SHORTWAVE
            IF (Aout(idSrad,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsrf(i,j,f)= FORCES(ng)%srflx(i,j)*fw(iw)
                END DO
              END DO
            END IF
#  endif
#  ifdef BULK_FLUXES
            IF (Aout(idLhea,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%fillhf(i,j,f)= FORCES(ng)%lhflx(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idShea,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%fillrf(i,j,f)= FORCES(ng)%shflx(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idLrad,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filshf(i,j,f)= FORCES(ng)%lrflx(i,j)*fw(iw)
                END DO
              END DO
            END IF
#   ifdef EMINUSP
            IF (Aout(idevap,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filevap(i,j,f)= FORCES(ng)%evap(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idrain,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filrain(i,j,f)= FORCES(ng)%rain(i,j)*fw(iw)
                END DO
              END DO
            END IF
#   endif
#  endif
# endif
            IF (Aout(idUsms,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsus(i,j,f) = FORCES(ng)%sustr(i,j)*     &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idVsms,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsvs(i,j,f) = FORCES(ng)%svstr(i,j)*     &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idUbms,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filbus(i,j,f) = FORCES(ng)%bustr(i,j)*     &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idVbms,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filbvs(i,j,f) = FORCES(ng)%bvstr(i,j)*     &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(id2dPV,ng)) THEN
              DO j=Jstr,JendR
                DO i=Istr,IendR
                  FILTER(ng)%filpotvor_bar(i,j,f) = potvor_bar(i,j)*    &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(id2dRV,ng)) THEN
              DO j=Jstr,JendR
                DO i=Istr,IendR
                  FILTER(ng)%filrelvor_bar(i,j,f) = relvor_bar(i,j)*    &
     &                          fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(id3dPV,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filpotvor(i,j,k,f) = potvor(i,j,k)*      &
     &                            fw(iw)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(id3dRV,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filrelvor(i,j,k,f) = relvor(i,j,k)*      &
     &                            fw(iw)
                  END DO
                END DO
              END DO
            END IF
# ifdef ICE_MODEL
            IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,Iend+1
                  FILTER(ng)%filui(i,j,f) = ICE(ng)%ui(i,j,Iuout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
              DO j=Jstr,Jend+1
                DO i=IstrR,IendR
                  FILTER(ng)%filvi(i,j,f) = ICE(ng)%vi(i,j,Iuout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idAice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filai(i,j,f) = ICE(ng)%ai(i,j,Iout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idHice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhi(i,j,f) = ICE(ng)%hi(i,j,Iout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idHsno,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhsn(i,j,f) = ICE(ng)%hsn(i,j,Iout)*     &
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idTice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filtis(i,j,f) = ICE(ng)%tis(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idTimid,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filti(i,j,f) = ICE(ng)%ti(i,j,Iout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idT2ice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filt2(i,j,f) = ICE(ng)%t2(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idAgeice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAgeice(i,j,f) =                         &
     &                            ICE(ng)%ageice(i,j,Iout)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idS0mk,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%fils0mk(i,j,f) = ICE(ng)%s0mk(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idT0mk,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filt0mk(i,j,f) = ICE(ng)%t0mk(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWfr,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWfr(i,j,f) = ICE(ng)%wfr(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWai,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWai(i,j,f) = ICE(ng)%wai(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWao,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWao(i,j,f) = ICE(ng)%wao(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWio,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWio(i,j,f) = ICE(ng)%wio(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWro,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWro(i,j,f) = ICE(ng)%wro(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idWdiv,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filWdiv(i,j,f) = ICE(ng)%wdiv(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSig11,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsig11(i,j,f) = ICE(ng)%sig11(i,j,Ieout)*&
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSig12,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsig12(i,j,f) = ICE(ng)%sig12(i,j,Ieout)*&
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSig22,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsig22(i,j,f) = ICE(ng)%sig22(i,j,Ieout)*&
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idIomflx,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filio_mflux(i,j,f) = ICE(ng)%io_mflux(i,j)*&
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idTauiw,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filutau_iw(i,j,f) = ICE(ng)%utau_iw(i,j)*  &
     &                            fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idChuiw,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filchu_iw(i,j,f) = ICE(ng)%chu_iw(i,j)*    &
     &                            fw(iw)
                END DO
              END DO
            END IF
#  ifdef ICE_DIAGS
            IF (Aout(idAlbe,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAlbe(i,j,f) =            &
     &                        FORCES(ng)%albedo(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idAlbI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAlbI(i,j,f) =            &
     &                        FORCES(ng)%albedo_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idLWdn,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filLWdn(i,j,f) =            &
     &                        FORCES(ng)%LW_down(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSWdn,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSWdn(i,j,f) =            &
     &                        FORCES(ng)%SW_down(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idLatI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filLatI(i,j,f) =            &
     &                        FORCES(ng)%lat_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSenI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSenI(i,j,f) =            &
     &                        FORCES(ng)%sens_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idLWupI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filLWupI(i,j,f) =            &
     &                        FORCES(ng)%LW_up_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSWupI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSWupI(i,j,f) =            &
     &                        FORCES(ng)%SW_up_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSWinI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSWinI(i,j,f) =            &
     &                        FORCES(ng)%qi_o_n(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSWthruI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSWthruI(i,j,f) =            &
     &                        FORCES(ng)%SW_thru_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSaltI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSaltI(i,j,f) =            &
     &                        FORCES(ng)%saltflux_ice(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSaltL,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSaltL(i,j,f) =            &
     &                        FORCES(ng)%saltflux_ocean(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idQao,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filQao(i,j,f) =            &
     &                        FORCES(ng)%qao_n(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idQai,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filQai(i,j,f) =            &
     &                        FORCES(ng)%qai_n(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idQio,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filQio(i,j,f) =            &
     &                        FORCES(ng)%qio_n(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idQi2,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filQi2(i,j,f) =            &
     &                        FORCES(ng)%qi2_n(i,j)*fw(iw)
                END DO
              END DO
            END IF
            IF (Aout(idSnoI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filSnoI(i,j,f) =            &
     &                        FORCES(ng)%snoice(i,j)*fw(iw)
                END DO
              END DO
            END IF
#  endif
# endif
            fcount(f)=iw
!
!-----------------------------------------------------------------------
!  Accumulate time-filtered fields.
!-----------------------------------------------------------------------
!
          ELSE
            iw=fcount(f)+1
             IF (Master) THEN
               print *, 'SET_FILTER add to sum', f, '    iw:', iw
             END IF
            IF (fcount(f).gt.0) THEN
              IF (Aout(idFsur,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filzeta(i,j,f)=FILTER(ng)%filzeta(i,j,f)+&
     &                               OCEAN(ng)%zeta(i,j,Kout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
                DO j=JstrR,JendR
                  DO i=Istr,Iend+1
                    FILTER(ng)%filu2d (i,j,f)=FILTER(ng)%filu2d (i,j,f)+&
     &                               OCEAN(ng)%ubar(i,j,Kout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
                DO j=Jstr,Jend+1
                  DO i=IstrR,IendR
                    FILTER(ng)%filv2d (i,j,f)=FILTER(ng)%filv2d (i,j,f)+&
     &                               OCEAN(ng)%vbar(i,j,Kout)*fw(iw)
                  END DO
                END DO
              END IF
# ifdef SOLVE3D
              IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
                DO k=1,N(ng)
                  DO j=JstrR,JendR
                    DO i=Istr,Iend+1
                      FILTER(ng)%filu3d(i,j,k,f)=                       &
     &                                    FILTER(ng)%filu3d(i,j,k,f)+   &
     &                                    OCEAN(ng)%u(i,j,k,Nout)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
                DO k=1,N(ng)
                  DO j=Jstr,Jend+1
                    DO i=IstrR,IendR
                      FILTER(ng)%filv3d(i,j,k,f)=                       &
     &                                    FILTER(ng)%filv3d(i,j,k,f)+   &
     &                                    OCEAN(ng)%v(i,j,k,Nout)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idDano,ng)) THEN
                DO k=1,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filrho(i,j,k,f)=                       &
     &                                    FILTER(ng)%filrho(i,j,k,f)+   &
     &                                       OCEAN(ng)%rho(i,j,k)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              DO itrc=1,NAT
                IF (Aout(idTvar(itrc),ng)) THEN
                  DO k=1,N(ng)
                    DO j=JstrR,JendR
                      DO i=IstrR,IendR
                        FILTER(ng)%filt(i,j,k,itrc,f)=                  &
     &                           FILTER(ng)%filt(i,j,k,itrc,f)+         &
     &                           OCEAN(ng)%t(i,j,k,Nout,itrc)*fw(iw)
                      END DO
                    END DO
                  END DO
                END IF
              END DO
              IF (Aout(idOvel,ng)) THEN
                DO k=0,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filomega(i,j,k,f)=                     &
     &                                  FILTER(ng)%filomega(i,j,k,f)+   &
     &                                  OCEAN(ng)%W(i,j,k)*fw(iw)*      &
     &                                  GRID(ng)%pm(i,j)*               &
     &                                  GRID(ng)%pn(i,j)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idWvel,ng)) THEN
                DO k=0,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filw3d(i,j,k,f)=                       &
     &                                  FILTER(ng)%filw3d(i,j,k,f)+     &
     &                                  OCEAN(ng)%wvel(i,j,k)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
#  ifdef LMD_SKPP
              IF (Aout(idHsbl,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filhsbl(i,j,f)=FILTER(ng)%filhsbl(i,j,f)+&
     &                               MIXING(ng)%hsbl(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#  endif
#  ifdef LMD_BKPP
              IF (Aout(idHbbl,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filhbbl(i,j,f)=FILTER(ng)%filhbbl(i,j,f)+&
     &                               MIXING(ng)%hbbl(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#  endif
              IF (Aout(idVvis,ng)) THEN
                DO k=0,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filAKv(i,j,k,f)=                       &
     &                                      FILTER(ng)%filAKv(i,j,k,f)+ &
     &                                      MIXING(ng)%Akv(i,j,k)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idTdif,ng)) THEN
                DO k=0,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filAKt(i,j,k,f)=                       &
     &                                FILTER(ng)%filAKt(i,j,k,f)+       &
     &                                MIXING(ng)%Akt(i,j,k,itemp)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idSdif,ng)) THEN
                DO k=0,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filAKs(i,j,k,f)=                       &
     &                                FILTER(ng)%filAKs(i,j,k,f)+       &
     &                                MIXING(ng)%Akt(i,j,k,isalt)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(idTsur(itemp),ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filstf(i,j,f)=FILTER(ng)%filstf(i,j,f)+  &
     &                                FORCES(ng)%stflx(i,j,itemp)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idTsur(isalt),ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filswf(i,j,f)=FILTER(ng)%filswf(i,j,f)+  &
     &                                FORCES(ng)%stflx(i,j,isalt)*fw(iw)
                  END DO
                END DO
              END IF
#  ifdef BULK_FLUXES
              IF (Aout(idLhea,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%fillhf(i,j,f)=FILTER(ng)%fillhf(i,j,f)+  &
     &                                      FORCES(ng)%lhflx(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idLrad,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%fillrf(i,j,f)=FILTER(ng)%fillrf(i,j,f)+  &
     &                                      FORCES(ng)%shflx(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idLhea,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filshf(i,j,f)=FILTER(ng)%filshf(i,j,f)+  &
     &                                      FORCES(ng)%lrflx(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#   ifdef EMINUSP
              IF (Aout(idevap,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filevap(i,j,f)=FILTER(ng)%filevap(i,j,f)+&
     &                                      FORCES(ng)%evap(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idrain,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filrain(i,j,f)=FILTER(ng)%filrain(i,j,f)+&
     &                                      FORCES(ng)%rain(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#   endif
#  endif
#  ifdef SHORTWAVE
              IF (Aout(idSrad,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filsrf(i,j,f)=FILTER(ng)%filsrf(i,j,f)+  &
     &                                      FORCES(ng)%srflx(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#  endif
# endif
              IF (Aout(idUsms,ng)) THEN
                DO j=JstrR,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filsus(i,j,f)=FILTER(ng)%filsus(i,j,f)+  &
     &                                      FORCES(ng)%sustr(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idVsms,ng)) THEN
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filsvs(i,j,f)=FILTER(ng)%filsvs(i,j,f)+  &
     &                                      FORCES(ng)%svstr(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idUbms,ng)) THEN
                DO j=JstrR,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filbus(i,j,f)=FILTER(ng)%filbus(i,j,f)+  &
     &                                      FORCES(ng)%bustr(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idVbms,ng)) THEN
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filbvs(i,j,f)=FILTER(ng)%filbvs(i,j,f)+  &
     &                                      FORCES(ng)%bvstr(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(id2dPV,ng)) THEN
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filpotvor_bar(i,j,f)=                    &
     &                         FILTER(ng)%filpotvor_bar(i,j,f)+         &
     &                         potvor_bar(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(id2dRV,ng)) THEN
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    FILTER(ng)%filrelvor_bar(i,j,f)=                    &
     &                         FILTER(ng)%filrelvor_bar(i,j,f)+         &
     &                         relvor_bar(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(id3dPV,ng)) THEN
                DO k=1,N(ng)
                  DO j=Jstr,JendR
                    DO i=Istr,IendR
                      FILTER(ng)%filpotvor(i,j,k,f)=                    &
     &                           FILTER(ng)%filpotvor(i,j,k,f)+         &
     &                           potvor(i,j,k)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
              IF (Aout(id3dRV,ng)) THEN
                DO k=1,N(ng)
                  DO j=Jstr,JendR
                    DO i=Istr,IendR
                      FILTER(ng)%filrelvor(i,j,k,f)=                    &
     &                           FILTER(ng)%filrelvor(i,j,k,f)+         &
     &                           relvor(i,j,k)*fw(iw)
                    END DO
                  END DO
                END DO
              END IF
# ifdef ICE_MODEL
              IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
                DO j=JstrR,JendR
                  DO i=Istr,Iend+1
                    FILTER(ng)%filui(i,j,f)=FILTER(ng)%filui(i,j,f)+    &
     &                                    ICE(ng)%ui(i,j,Iuout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
                DO j=Jstr,Jend+1
                  DO i=IstrR,IendR
                    FILTER(ng)%filvi(i,j,f)=FILTER(ng)%filvi(i,j,f)+    &
     &                                    ICE(ng)%vi(i,j,Iuout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idAice,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filai(i,j,f)=FILTER(ng)%filai(i,j,f)+    &
     &                                    ICE(ng)%ai(i,j,Iout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idHice,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filhi(i,j,f)=FILTER(ng)%filhi(i,j,f)+    &
     &                                    ICE(ng)%hi(i,j,Iout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idHsno,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filhsn(i,j,f)=FILTER(ng)%filhsn(i,j,f)+  &
     &                                    ICE(ng)%hsn(i,j,Iout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idTice,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filtis(i,j,f)=FILTER(ng)%filtis(i,j,f)+  &
     &                                    ICE(ng)%tis(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idAgeice,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAgeice(i,j,f)=                        &
     &                              FILTER(ng)%filAgeice(i,j,f)+        &
     &                              ICE(ng)%ageice(i,j,Iout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idTimid,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filti(i,j,f)=FILTER(ng)%filti(i,j,f)+    &
     &                                    ICE(ng)%ti(i,j,Iout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idT2ice,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filt2(i,j,f)=FILTER(ng)%filt2(i,j,f)+    &
     &                                    ICE(ng)%t2(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idS0mk,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%fils0mk(i,j,f)=FILTER(ng)%fils0mk(i,j,f)+&
     &                                    ICE(ng)%s0mk(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idT0mk,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filt0mk(i,j,f)=FILTER(ng)%filt0mk(i,j,f)+&
     &                                    ICE(ng)%t0mk(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWfr,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWfr(i,j,f)=FILTER(ng)%filWfr(i,j,f)+  &
     &                                    ICE(ng)%wfr(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWai,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWai(i,j,f)=FILTER(ng)%filWai(i,j,f)+  &
     &                                    ICE(ng)%wai(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWao,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWao(i,j,f)=FILTER(ng)%filWao(i,j,f)+  &
     &                                    ICE(ng)%wao(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWio,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWio(i,j,f)=FILTER(ng)%filWio(i,j,f)+  &
     &                                    ICE(ng)%wio(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWro,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWro(i,j,f)=FILTER(ng)%filWro(i,j,f)+  &
     &                                    ICE(ng)%wro(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idWdiv,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filWdiv(i,j,f)=FILTER(ng)%filWdiv(i,j,f)+&
     &                                    ICE(ng)%wdiv(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSig11,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filsig11(i,j,f)=                         &
     &                               FILTER(ng)%filsig11(i,j,f)+        &
     &                               ICE(ng)%sig11(i,j,Ieout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSig12,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filsig12(i,j,f)=                         &
     &                               FILTER(ng)%filsig12(i,j,f)+        &
     &                               ICE(ng)%sig12(i,j,Ieout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSig22,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filsig22(i,j,f)=                         &
     &                               FILTER(ng)%filsig22(i,j,f)+        &
     &                               ICE(ng)%sig22(i,j,Ieout)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idIomflx,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filio_mflux(i,j,f)=                      &
     &                                FILTER(ng)%filio_mflux(i,j,f)+    &
     &                                    ICE(ng)%io_mflux(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idTauiw,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filutau_iw(i,j,f)=                       &
     &                                FILTER(ng)%filutau_iw(i,j,f)+     &
     &                                    ICE(ng)%utau_iw(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idChuiw,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filchu_iw(i,j,f)=                        &
     &                                FILTER(ng)%filchu_iw(i,j,f)+      &
     &                                    ICE(ng)%chu_iw(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#   ifdef ICE_DIAGS
              IF (Aout(idAlbe,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAlbe(i,j,f)=FILTER(ng)%filAlbe(i,j,f)+&
     &                                    FORCES(ng)%albedo(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idAlbI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAlbI(i,j,f)=FILTER(ng)%filAlbI(i,j,f)+&
     &                            FORCES(ng)%albedo_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idLWdn,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filLWdn(i,j,f)=FILTER(ng)%filLWdn(i,j,f)+&
     &                                    FORCES(ng)%LW_down(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSWdn,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSWdn(i,j,f)=FILTER(ng)%filSWdn(i,j,f)+&
     &                                    FORCES(ng)%SW_down(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idLatI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filLatI(i,j,f)=FILTER(ng)%filLatI(i,j,f)+&
     &                                    FORCES(ng)%lat_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSenI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSenI(i,j,f)=FILTER(ng)%filSenI(i,j,f)+&
     &                                   FORCES(ng)%sens_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idLWupI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filLWupI(i,j,f)=                         &
     &                                       FILTER(ng)%filLWupI(i,j,f)+&
     &                                  FORCES(ng)%LW_up_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSWupI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSWupI(i,j,f)=                         &
     &                                       FILTER(ng)%filSWupI(i,j,f)+&
     &                                 FORCES(ng)%SW_up_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSWinI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSWinI(i,j,f)=                         &
     &                                       FILTER(ng)%filSWinI(i,j,f)+&
     &                                FORCES(ng)%qi_o_n(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSWthruI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSWthruI(i,j,f)=                       &
     &                                     FILTER(ng)%filSWthruI(i,j,f)+&
     &                               FORCES(ng)%SW_thru_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSaltI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSaltI(i,j,f)=                         &
     &                                       FILTER(ng)%filSaltI(i,j,f)+&
     &                            FORCES(ng)%saltflux_ice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSaltL,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSaltL(i,j,f)=                         &
     &                                       FILTER(ng)%filSaltL(i,j,f)+&
     &                           FORCES(ng)%saltflux_ocean(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idQao,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filQao(i,j,f)=FILTER(ng)%filQao(i,j,f)+  &
     &                                    FORCES(ng)%qao_n(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idQai,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filQai(i,j,f)=FILTER(ng)%filQai(i,j,f)+  &
     &                                    FORCES(ng)%qai_n(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idQio,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filQio(i,j,f)=FILTER(ng)%filQio(i,j,f)+  &
     &                                    FORCES(ng)%qio_n(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idQi2,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filQi2(i,j,f)=FILTER(ng)%filQi2(i,j,f)+  &
     &                                    FORCES(ng)%qi2_n(i,j)*fw(iw)
                  END DO
                END DO
              END IF
              IF (Aout(idSnoI,ng)) THEN
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filSnoI(i,j,f)=FILTER(ng)%filSnoI(i,j,f)+&
     &                                    FORCES(ng)%snoice(i,j)*fw(iw)
                  END DO
                END DO
              END IF
#   endif
# endif
            END IF
            fcount(f)=iw
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Transfer finished sums into averages arrays for writing out.
!-----------------------------------------------------------------------
!
      IF ((MOD(iic(ng)-1,navg(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)))) THEN
!    &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        DO f=1,FILN
          IF (fcount(f).eq.nw) THEN
            IF (Master) THEN
              print *,'SET_FILTER transferring sums', f
            END IF
            IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
              avgtime(ng)=time(ng)-0.5_r8*REAL(nfil,r8)*dt(ng)
            END IF
            IF (Aout(idFsur,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgzeta(i,j)=FILTER(ng)%filzeta(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idUbar,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,Iend+1
                  AVERAGE(ng)%avgu2d (i,j)=FILTER(ng)%filu2d (i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idVbar,ng)) THEN
              DO j=Jstr,Jend+1
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgv2d (i,j)=FILTER(ng)%filv2d (i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
              CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      FILTER(ng) % filu2d(:,:,f),                 &
     &                      FILTER(ng) % filv2d(:,:,f),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
            END IF
# ifdef SOLVE3D
            IF (Aout(idUvel,ng)) THEN
              DO k=1,N(ng)
                DO j=JstrR,JendR
                  DO i=Istr,Iend+1
                    AVERAGE(ng)%avgu3d(i,j,k)=FILTER(ng)%filu3d(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idVvel,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgv3d(i,j,k)=FILTER(ng)%filv3d(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
              CALL uv_rotate3d (ng, tile, .FALSE., .FALSE.,             &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      FILTER(ng) % filu3d(:,:,:,f),               &
     &                      FILTER(ng) % filv3d(:,:,:,f),               &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
            END IF
            IF (Aout(idDano,ng)) THEN
              DO k=1,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgrho(i,j,k)=FILTER(ng)%filrho(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            DO itrc=1,NAT
              IF (Aout(idTvar(itrc),ng)) THEN
                DO k=1,N(ng)
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      AVERAGE(ng)%avgt(i,j,k,itrc)=                     &
     &                                FILTER(ng)%filt(i,j,k,itrc,f)
                    END DO
                  END DO
                END DO
              END IF
            END DO
            IF (Aout(idOvel,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgw3d(i,j,k)=                          &
     &                          FILTER(ng)%filomega(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idWvel,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgwvel(i,j,k)=                         &
     &                          FILTER(ng)%filw3d(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
#  ifdef LMD_SKPP
            IF (Aout(idHsbl,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avghsbl(i,j)=FILTER(ng)%filhsbl(i,j,f)
                END DO
              END DO
            END IF
#  endif
#  ifdef LMD_BKPP
            IF (Aout(idHbbl,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avghbbl(i,j)=FILTER(ng)%filhbbl(i,j,f)
                END DO
              END DO
            END IF
#  endif
            IF (Aout(idVvis,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgAKv(i,j,k)=FILTER(ng)%filAKv(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idTdif,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgAKt(i,j,k)=FILTER(ng)%filAKt(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idSdif,ng)) THEN
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgAKs(i,j,k)=FILTER(ng)%filAKs(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(idTsur(itemp),ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgstf(i,j)=FILTER(ng)%filstf(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idTsur(isalt),ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgswf(i,j)=FILTER(ng)%filswf(i,j,f)
                END DO
              END DO
            END IF
#  ifdef BULK_FLUXES
            IF (Aout(idLhea,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avglhf(i,j)=FILTER(ng)%fillhf(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idLrad,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avglrf(i,j)=FILTER(ng)%fillrf(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idShea,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgshf(i,j)=FILTER(ng)%filshf(i,j,f)
                END DO
              END DO
            END IF
#   ifdef EMINUSP
            IF (Aout(idevap,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgevap(i,j)=FILTER(ng)%filevap(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idrain,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgrain(i,j)=FILTER(ng)%filrain(i,j,f)
                END DO
              END DO
            END IF
#   endif
#  endif
#  ifdef SHORTWAVE
            IF (Aout(idSrad,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgsrf(i,j)=FILTER(ng)%filsrf(i,j,f)
                END DO
              END DO
            END IF
#  endif
# endif
            IF (Aout(idUsms,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  AVERAGE(ng)%avgsus(i,j)=FILTER(ng)%filsus(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idVsms,ng)) THEN
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgsvs(i,j)=FILTER(ng)%filsvs(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idUbms,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  AVERAGE(ng)%avgbus(i,j)=FILTER(ng)%filbus(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idVbms,ng)) THEN
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgbvs(i,j)=FILTER(ng)%filbvs(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(id2dPV,ng)) THEN
              DO j=Jstr,JendR
                DO i=Istr,IendR
                  AVERAGE(ng)%avgpvor2d(i,j)=                           &
     &                        FILTER(ng)%filpotvor_bar(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(id2dRV,ng)) THEN
              DO j=Jstr,JendR
                DO i=Istr,IendR
                  AVERAGE(ng)%avgrvor2d(i,j)=                           &
     &                        FILTER(ng)%filrelvor_bar(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(id3dPV,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    AVERAGE(ng)%avgpvor3d(i,j,k)=                       &
     &                          FILTER(ng)%filpotvor(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
            IF (Aout(id3dRV,ng)) THEN
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=Istr,IendR
                    AVERAGE(ng)%avgrvor3d(i,j,k)=                       &
     &                          FILTER(ng)%filrelvor(i,j,k,f)
                  END DO
                END DO
              END DO
            END IF
# ifdef ICE_MODEL
            IF (Aout(idUice,ng)) THEN
              DO j=JstrR,JendR
                DO i=Istr,Iend+1
                  AVERAGE(ng)%avguice(i,j)=FILTER(ng)%filui(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idVice,ng)) THEN
              DO j=Jstr,Jend+1
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgvice(i,j)=FILTER(ng)%filvi(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
              CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      FILTER(ng) % filui(:,:,f),                  &
     &                      FILTER(ng) % filvi(:,:,f),                  &
     &                      AVERAGE(ng)%avguiceE,                       &
     &                      AVERAGE(ng)%avgviceN)
            END IF
            IF (Aout(idAice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgaice(i,j)=FILTER(ng)%filai(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idHice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avghice(i,j)=FILTER(ng)%filhi(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idHsno,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avghsno(i,j)=FILTER(ng)%filhsn(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idTice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgtice(i,j)=FILTER(ng)%filtis(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idTimid,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgtimid(i,j)=FILTER(ng)%filti(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idT2ice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt2ice(i,j)=FILTER(ng)%filt2(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idAgeice,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAgeice(i,j)=FILTER(ng)%filAgeice(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idS0mk,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgS0mk(i,j)=FILTER(ng)%fils0mk(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idT0mk,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgT0mk(i,j)=FILTER(ng)%filt0mk(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWfr,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWfr(i,j)=FILTER(ng)%filWfr(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWai,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWai(i,j)=FILTER(ng)%filWai(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWao,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWao(i,j)=FILTER(ng)%filWao(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWio,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWio(i,j)=FILTER(ng)%filWio(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWro,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWro(i,j)=FILTER(ng)%filWro(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idWdiv,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgWdiv(i,j)=FILTER(ng)%filWdiv(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSig11,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgsig11(i,j)=FILTER(ng)%filsig11(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSig12,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgsig12(i,j)=FILTER(ng)%filsig12(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSig22,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgsig22(i,j)=FILTER(ng)%filsig22(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idIomflx,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgiomflx(i,j)=                           &
     &                          FILTER(ng)%filio_mflux(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idTauiw,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgutau_iw(i,j)=                          &
     &                          FILTER(ng)%filutau_iw(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idChuiw,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgchu_iw(i,j)=FILTER(ng)%filchu_iw(i,j,f)
                END DO
              END DO
            END IF
#   ifdef ICE_DIAGS
            IF (Aout(idAlbe,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAlbe(i,j)=FILTER(ng)%filAlbe(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idAlbI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAlbI(i,j)=FILTER(ng)%filAlbI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idLWdn,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgLWdn(i,j)=FILTER(ng)%filLWdn(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSWdn,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSWdn(i,j)=FILTER(ng)%filSWdn(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idLatI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgLatI(i,j)=FILTER(ng)%filLatI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSenI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSenI(i,j)=FILTER(ng)%filSenI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idLWupI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgLWupI(i,j)=FILTER(ng)%filLWupI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSWupI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSWupI(i,j)=FILTER(ng)%filSWupI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSWinI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSWinI(i,j)=FILTER(ng)%filSWinI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSWthruI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSWthruI(i,j)=                          &
     &                          FILTER(ng)%filSWthruI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSaltI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSaltI(i,j)=FILTER(ng)%filSaltI(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSaltL,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSaltL(i,j)=FILTER(ng)%filSaltL(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idQao,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgQao(i,j)=FILTER(ng)%filQao(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idQai,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgQai(i,j)=FILTER(ng)%filQai(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idQio,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgQio(i,j)=FILTER(ng)%filQio(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idQi2,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgQi2(i,j)=FILTER(ng)%filQi2(i,j,f)
                END DO
              END DO
            END IF
            IF (Aout(idSnoI,ng)) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgSnoI(i,j)=FILTER(ng)%filSnoI(i,j,f)
                END DO
              END DO
            END IF
#   endif
# endif
            fcount(f)=0
          END IF
        END DO
      END IF
!     IF (Master) print *, 'end of SET_FILTER'
      RETURN
      END SUBROUTINE set_filter_tile
#else
      SUBROUTINE set_filter
      RETURN
      END SUBROUTINE set_filter
#endif
