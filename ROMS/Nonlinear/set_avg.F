#include "cppdefs.h"
      MODULE set_avg_mod
#ifdef AVERAGES
!
!svn $Id: set_avg.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!  It computes leas-squares coefficients to detide time-averaged       !
!  fields. Notice that "set_detide" is called last since we need       !
!  the regular time-averages for those fields to detide.               !
!                                                                      !
# endif
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef WET_DRY
      USE mod_grid
# endif
      USE mod_stepping
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
# ifdef WET_DRY
!
      USE set_masks_mod, ONLY : set_avg_masks
# endif
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
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
#  endif
     &                   NOUT,                                          &
# endif
     &                   KOUT)

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      CALL set_detide_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NTC(ng), KOUT,                              &
#  ifdef SOLVE3D
     &                      NOUT,                                       &
#  endif
     &                      TIDES(ng) % CosOmega,                       &
     &                      TIDES(ng) % SinOmega,                       &
     &                      TIDES(ng) % CosW_avg,                       &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      TIDES(ng) % SinW_avg,                       &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      TIDES(ng) % SinWCosW)
# endif

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
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
#  endif
     &                         Nout,                                    &
# endif
     &                         Kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
# ifdef SOLVE3D
      USE mod_mixing
#  ifdef ICE_MODEL
      USE mod_ice
#  endif
# endif
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT && defined BEDLOAD
      USE mod_sedbed
      USE mod_sediment
# endif
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
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
      integer, intent(in) :: Kout
# ifdef SOLVE3D
#  ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
#  endif
      integer, intent(in) :: Nout
# endif
!
!  Local variable declarations.
!
      integer :: i, it, itrc, j, k
# ifdef  BIOFLUX
       integer ::itrc2
# endif

      real(r8) :: fac

      real(r8) :: pfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

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
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity fields.
!-----------------------------------------------------------------------
!
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
!
!  Initialize state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ubar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%vbar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)*    &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)*    &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
#  ifdef WET_DRY
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=OCEAN(ng)%wvel(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)*  &
     &                                     GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)*&
     &                                       GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO

#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=OCEAN(ng)%PONsed(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgPONsed(i,j)=AVERAGE(ng)%avgPONsed(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=OCEAN(ng)%OPALsed(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALsed(i,j)*  &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=OCEAN(ng)%DENITsed(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgDENITsed(i,j)=AVERAGE(ng)%avgDENITsed(i,j)*&
     &                                     GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=OCEAN(ng)%PON_burial(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgPONbur(i,j)=AVERAGE(ng)%avgPONbur(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALbur(i,j)=OCEAN(ng)%OPAL_burial(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgOPALbur(i,j)=AVERAGE(ng)%avgOPALbur(i,j)*  &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (Aout(idNPP,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgNPP(i,j)=OCEAN(ng)%Bio_NPP(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgNPP(i,j)=AVERAGE(ng)%avgNPP(i,j)*  &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%bedldu(i,j,it)
#   ifdef WET_DRY
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%bedldu(i,j,it)* &
     &                                       GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%bedldv(i,j,it)
#   ifdef WET_DRY
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%bedldv(i,j,it)* &
     &                                       GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END IF
        END DO
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%Akt(i,j,k,isalt)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Initialize 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=COUPLING(ng)%DU_avg1(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)*  &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=COUPLING(ng)%DU_avg2(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)*  &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=COUPLING(ng)%DV_avg1(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)*  &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=COUPLING(ng)%DV_avg2(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)*  &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF
# ifdef BBL
        IF (Aout(idUbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbrs(i,j)=FORCES(ng)%bustrc(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUbrs(i,j)=AVERAGE(ng)%avgUbrs(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbrs(i,j)=FORCES(ng)%bvstrc(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVbrs(i,j)=AVERAGE(ng)%avgVbrs(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbws(i,j)=FORCES(ng)%bustrw(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUbws(i,j)=AVERAGE(ng)%avgUbws(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbws(i,j)=FORCES(ng)%bvstrw(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVbws(i,j)=AVERAGE(ng)%avgVbws(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbcs(i,j)=FORCES(ng)%bustrcwmax(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUbcs(i,j)=AVERAGE(ng)%avgUbcs(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbcs(i,j)=FORCES(ng)%bvstrcwmax(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVbcs(i,j)=AVERAGE(ng)%avgVbcs(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUVwc,ng)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj))
          wrk(LBi:UBi,LBj:UBj)=0.0_r8
          scale=rho0
          gtype=gfactor*r2dvar
          wrk=sqrt(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+               &
     &             BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUVwc(i,j)=wrk(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUVwc(i,j)=AVERAGE(ng)%avgUVwc(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
          deallocate(wrk)
        END IF
        IF (Aout(idUbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbot(i,j)=BBL(ng)%Ubot(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUbot(i,j)=AVERAGE(ng)%avgUbot(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbot(i,j)=BBL(ng)%Vbot(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVbot(i,j)=AVERAGE(ng)%avgVbot(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbur(i,j)=BBL(ng)%Ubur(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUbur(i,j)=AVERAGE(ng)%avgUbur(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVbvr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbvr(i,j)=BBL(ng)%Vbur(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVbvr(i,j)=AVERAGE(ng)%avgVbvr(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef SOLVE3D
#  if defined RUNOFF
        IF (Aout(idRunoff,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgRunoff(i,j)=FORCES(ng)%Runoff(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgRunoff(i,j)=AVERAGE(ng)%avgRunoff(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=FORCES(ng)%Pair(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgPair(i,j)=AVERAGE(ng)%avgPair(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=FORCES(ng)%Uwind(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgUwind(i,j)=AVERAGE(ng)%avgUwind(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=FORCES(ng)%Vwind(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgVwind(i,j)=AVERAGE(ng)%avgVwind(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#   ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#   endif
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      AVERAGE(ng)%avgUwindE,                      &
     &                      AVERAGE(ng)%avgVwindN)
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
#  ifdef WET_DRY
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
#  ifdef WET_DRY
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=FORCES(ng)%sssflx(i,j,isalt)
            END DO
          END DO
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
#    ifdef WET_DRY
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
#    ifdef WET_DRY
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=ICE(ng)%ui(i,j,Iuout)
#  ifdef WET_DRY
              AVERAGE(ng)%avguice(i,j)=AVERAGE(ng)%avguice(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=ICE(ng)%vi(i,j,Iuout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgvice(i,j)=AVERAGE(ng)%avgvice(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      ICE(ng) % ui(:,:,Kout),                     &
     &                      ICE(ng) % vi(:,:,Kout),                     &
     &                      AVERAGE(ng)%avguiceE,                       &
     &                      AVERAGE(ng)%avgviceN)
        END IF

        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=ICE(ng)%ai(i,j,Iout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgaice(i,j)=AVERAGE(ng)%avgaice(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=ICE(ng)%hi(i,j,Iout)
#  ifdef WET_DRY
              AVERAGE(ng)%avghice(i,j)=AVERAGE(ng)%avghice(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=ICE(ng)%hsn(i,j,Iout)
#  ifdef WET_DRY
              AVERAGE(ng)%avghsno(i,j)=AVERAGE(ng)%avghsno(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=ICE(ng)%tis(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgtice(i,j)=AVERAGE(ng)%avgtice(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=ICE(ng)%ti(i,j,Iout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgtimid(i,j)=AVERAGE(ng)%avgtimid(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
#  ifdef MELT_PONDS
        IF (Aout(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgapond(i,j)=ICE(ng)%apond(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avgapond(i,j)=AVERAGE(ng)%avgapond(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghpond(i,j)=ICE(ng)%hpond(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avghpond(i,j)=AVERAGE(ng)%avghpond(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=ICE(ng)%ageice(i,j,Iout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgageice(i,j)=AVERAGE(ng)%avgageice(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=ICE(ng)%io_mflux(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgiomflx(i,j)=AVERAGE(ng)%avgiomflx(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=ICE(ng)%sig11(i,j,Ieout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgsig11(i,j)=AVERAGE(ng)%avgsig11(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=ICE(ng)%sig12(i,j,Ieout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgsig12(i,j)=AVERAGE(ng)%avgsig12(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=ICE(ng)%sig22(i,j,Ieout)
#  ifdef WET_DRY
              AVERAGE(ng)%avgsig22(i,j)=AVERAGE(ng)%avgsig22(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=ICE(ng)%t0mk(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgT0mk(i,j)=AVERAGE(ng)%avgT0mk(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=ICE(ng)%s0mk(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgS0mk(i,j)=AVERAGE(ng)%avgS0mk(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=ICE(ng)%wfr(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWfr(i,j)=AVERAGE(ng)%avgWfr(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=ICE(ng)%wai(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWai(i,j)=AVERAGE(ng)%avgWai(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=ICE(ng)%wao(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWao(i,j)=AVERAGE(ng)%avgWao(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=ICE(ng)%wio(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWio(i,j)=AVERAGE(ng)%avgWio(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=ICE(ng)%wro(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWro(i,j)=AVERAGE(ng)%avgWro(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWdiv,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiv(i,j)=ICE(ng)%wdiv(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdiv(i,j)=AVERAGE(ng)%avgWdiv(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=ICE(ng)%utau_iw(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgutau_iw(i,j)=AVERAGE(ng)%avgutau_iw(i,j)*  &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=ICE(ng)%chu_iw(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgchu_iw(i,j)=AVERAGE(ng)%avgchu_iw(i,j)*    &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef ICE_BIO
#  ifdef BERING_10K
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINO3(i,j)=ICE(ng)%IceNO3(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avgINO3(i,j)=AVERAGE(ng)%avgINO3(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgIPhL(i,j)=ICE(ng)%IcePhL(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avgIPhL(i,j)=AVERAGE(ng)%avgIPhL(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINH4(i,j)=ICE(ng)%IceNH4(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avgINH4(i,j)=AVERAGE(ng)%avgINH4(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgILog(i,j)=ICE(ng)%IceLog(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE(ng)%avgILog(i,j)=AVERAGE(ng)%avgILog(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BEST_NPZ
#  ifdef STATIONARY
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                       OCEAN(ng)%st(i,j,k,Nout,itrc)
#   ifdef WET_DRY
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                  AVERAGE(ng)%avgst(i,j,k,itrc)*                  &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef PROD3
        DO itrc=1,NPT3(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgpt3(i,j,k,itrc)=                         &
     &                      OCEAN(ng)%pt3(i,j,k,Nout,itrc)
#   ifdef WET_DRY
                AVERAGE(ng)%avgpt3(i,j,k,itrc)=                         &
     &                  AVERAGE(ng)%avgpt3(i,j,k,itrc)*                 &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef PROD2
        DO itrc=1,NPT2(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgpt2(i,j,itrc)=                             &
     &                      OCEAN(ng)%pt2(i,j,Nout,itrc)
#   ifdef WET_DRY
              AVERAGE(ng)%avgpt2(i,j,itrc)=                             &
     &                  AVERAGE(ng)%avgpt2(i,j,itrc)*                   &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END DO
#  endif
#  ifdef STATIONARY2
        DO itrc=1,NTS2(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgst2(i,j,itrc)=                             &
     &                      OCEAN(ng)%st2(i,j,Nout,itrc)
#   ifdef WET_DRY
              AVERAGE(ng)%avgst2(i,j,itrc)=AVERAGE(ng)%avgst2(i,j,itrc)*&
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END DO
#  endif
#  ifdef BENTHIC
        DO itrc=1,NBeT(ng)
          DO k=1,NBL(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgbt(i,j,k,itrc)=                          &
     &                      OCEAN(ng)%bt(i,j,k,Nout,itrc)
#   ifdef WET_DRY
                AVERAGE(ng)%avgbt(i,j,k,itrc)=                          &
     &                  AVERAGE(ng)%avgpt(i,j,k,itrc)*                  &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END DO
#  endif
#  if defined CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgit(i,j,itrc)=                              &
     &                      OCEAN(ng)%it(i,j,Nout,itrc)
#   ifdef WET_DRY
              AVERAGE(ng)%avgit(i,j,itrc)=AVERAGE(ng)%avgit(i,j,itrc)*  &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END DO
#  endif
#  ifdef BIOFLUX
        DO itrc=1,NT(ng)
          DO itrc2=1,NT(ng)
            AVERAGE(ng)%avgbflx(itrc,itrc2)=                            &
     &                      OCEAN(ng)%bflx(itrc,itrc2)
          END DO
        END DO
#  endif
# endif

# ifdef WEC
!
!  Initialize WEC values.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=OCEAN(ng)%ubar_stokes(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=OCEAN(ng)%vbar_stokes(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=MIXING(ng)%rustr2d(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgu2rs(i,j)=AVERAGE(ng)%avgu2rs(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=MIXING(ng)%rvstr2d(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgv2rs(i,j)=AVERAGE(ng)%avgv2rs(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WEC_MELLOR
        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=MIXING(ng)%Sxx_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSxx2d(i,j)=AVERAGE(ng)%avgSxx2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=MIXING(ng)%Sxy_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSxy2d(i,j)=AVERAGE(ng)%avgSxy2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=MIXING(ng)%Syy_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSyy2d(i,j)=AVERAGE(ng)%avgSyy2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef SOLVE3D
#  ifdef WEC
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=OCEAN(ng)%u_stokes(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=OCEAN(ng)%v_stokes(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3Sd(i,j,k)=OCEAN(ng)%W_stokes(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgw3Sd(i,j,k)=AVERAGE(ng)%avgw3Sd(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3St(i,j,k)=OCEAN(ng)%wstvel(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgw3St(i,j,k)=AVERAGE(ng)%avgw3St(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=MIXING(ng)%rustr3d(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgu3rs(i,j,k)=AVERAGE(ng)%avgu3rs(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3rs(i,j,k)=MIXING(ng)%rvstr3d(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgv3rs(i,j,k)=AVERAGE(ng)%avgv3rs(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef WEC_VF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=OCEAN(ng)%zetaw(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgWztw(i,j)=AVERAGE(ng)%avgWztw(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgwqsp(i,j)=OCEAN(ng)%qsp(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgwqsp(i,j)=AVERAGE(ng)%avgwqsp(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgwbeh(i,j)=OCEAN(ng)%bh(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgwbeh(i,j)=AVERAGE(ng)%avgwbeh(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef WEC_MELLOR
        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=MIXING(ng)%Sxx(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSxx3d(i,j,k)=AVERAGE(ng)%avgSxx3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=MIXING(ng)%Sxy(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSxy3d(i,j,k)=AVERAGE(ng)%avgSxy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=MIXING(ng)%Syy(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSyy3d(i,j,k)=AVERAGE(ng)%avgSyy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=MIXING(ng)%Szx(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSzx3d(i,j,k)=AVERAGE(ng)%avgSzx3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=MIXING(ng)%Szy(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSzy3d(i,j,k)=AVERAGE(ng)%avgSzy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
# ifdef WAVES_HEIGHT
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=FORCES(ng)%Hwave(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWamp(i,j)=AVERAGE(ng)%avgWamp(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=FORCES(ng)%Lwave(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWlen(i,j)=AVERAGE(ng)%avgWlen(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (Aout(idWlep,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlep(i,j)=FORCES(ng)%Lwavep(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWlep(i,j)=AVERAGE(ng)%avgWlep(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_DIR
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=FORCES(ng)%Dwave(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdir(i,j)=AVERAGE(ng)%avgWdir(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=FORCES(ng)%Pwave_top(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWptp(i,j)=AVERAGE(ng)%avgWptp(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=FORCES(ng)%Pwave_bot(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWpbt(i,j)=AVERAGE(ng)%avgWpbt(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef BBL_MODEL
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=FORCES(ng)%Uwave_rms(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWorb(i,j)=AVERAGE(ng)%avgWorb(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=FORCES(ng)%Dissip_fric(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdif(i,j)=AVERAGE(ng)%avgWdif(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=FORCES(ng)%Dissip_break(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdib(i,j)=AVERAGE(ng)%avgWdib(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=FORCES(ng)%Dissip_wcap(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdiw(i,j)=AVERAGE(ng)%avgWdiw(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef ROLLER_SVENDSEN
        IF (Aout(idWbrk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbrk(i,j)=FORCES(ng)%Wave_break(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWbrk(i,j)=AVERAGE(ng)%avgWbrk(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef WEC_ROLLER
        IF (Aout(idWdis,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdis(i,j)=FORCES(ng)%Dissip_roller(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWdis(i,j)=AVERAGE(ng)%avgWdis(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idWrol,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWrol(i,j)=FORCES(ng)%rollA(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgWrol(i,j)=AVERAGE(ng)%avgWrol(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
# ifdef UV_KIRBY
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=OCEAN(ng)%uWave(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgUwav(i,j)=AVERAGE(ng)%avgUwav(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=OCEAN(ng)%vWave(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgVwav(i,j)=AVERAGE(ng)%avgVwav(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
!
!  Initialize vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=potvor_bar(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)*    &
     &                                   GRID(ng)%pmask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=relvor_bar(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)*    &
     &                                   GRID(ng)%pmask_full(i,j)
# endif
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=potvor(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)*  &
     &                                       GRID(ng)%pmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=relvor(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)*  &
     &                                       GRID(ng)%pmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Initialize quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)*            &
     &                               GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)*            &
     &                               GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)*            &
     &                               GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)*      &
     &                                   GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)*      &
     &                                   GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
#  ifdef WET_DRY
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)*      &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=GRID(ng)%Huon(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=GRID(ng)%Hvom(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%umask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%vmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT   &
     &                                                      (i,j,k,it)* &
     &                                          GRID(ng)%umask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT   &
     &                                                      (i,j,k,it)* &
     &                                          GRID(ng)%vmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
!
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
!
!  Accumulate state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%rmask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
#  ifdef WET_DRY
     &                                       GRID(ng)%rmask_full(i,j)*  &
#  endif
     &                                       OCEAN(ng)%t(i,j,k,Nout,it)
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
#   ifdef WET_DRY
     &                                   GRID(ng)%rmask_full(i,j)*      &
#   endif
     &                                   OCEAN(ng)%PONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALsed(i,j)+  &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 OCEAN(ng)%OPALsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=AVERAGE(ng)%avgDENITsed(i,j)+&
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 OCEAN(ng)%DENITsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=AVERAGE(ng)%avgPONbur(i,j)+    &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 OCEAN(ng)%PON_burial(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALbur(i,j)=AVERAGE(ng)%avgOPALbur(i,j)+  &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 OCEAN(ng)%OPAL_burial(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (Aout(idNPP,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgNPP(i,j)=AVERAGE(ng)%avgNPP(i,j)+          &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 OCEAN(ng)%Bio_NPP(i,j)
            END DO
          END DO
        END IF
#  endif
#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)+  &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       SEDBED(ng)%bedldu(i,j,it)
              END DO
            END DO
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)+  &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       SEDBED(ng)%bedldv(i,j,it)
              END DO
            END DO
          END IF
        END DO
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Accumulate 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
# ifdef BBL
        IF (Aout(idUbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbrs(i,j)=AVERAGE(ng)%avgUbrs(i,j)+        &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 FORCES(ng)%bustrc(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbrs(i,j)=AVERAGE(ng)%avgVbrs(i,j)+        &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 FORCES(ng)%bvstrc(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbws(i,j)=AVERAGE(ng)%avgUbws(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bustrw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbws(i,j)=AVERAGE(ng)%avgVbws(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bvstrw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbcs(i,j)=AVERAGE(ng)%avgUbcs(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bustrcwmax(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUVwc,ng)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj))
          wrk(LBi:UBi,LBj:UBj)=0.0_r8
          scale=rho0
          gtype=gfactor*r2dvar
          wrk=sqrt(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+               &
     &             BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUVwc(i,j)=AVERAGE(ng)%avgUVwc(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                wrk(i,j)
            END DO
          END DO
          deallocate(wrk)
        END IF
        IF (Aout(idUbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbot(i,j)=AVERAGE(ng)%avgUbot(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                BBL(ng)%Ubot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbot(i,j)=AVERAGE(ng)%avgVbot(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                BBL(ng)%Vbot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbur(i,j)=AVERAGE(ng)%avgUbur(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                BBL(ng)%Ubur(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbvr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbvr(i,j)=AVERAGE(ng)%avgVbvr(i,j)+        &
# ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
# endif
     &                                BBL(ng)%Vbur(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef SOLVE3D

#  if defined RUNOFF
        IF (Aout(idRunoff,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgRunoff(i,j)=AVERAGE(ng)%avgRunoff(i,j)+    &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 FORCES(ng)%Runoff(i,j)
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=AVERAGE(ng)%avgPair(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 FORCES(ng)%Pair(i,j)
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=AVERAGE(ng)%avgUwind(i,j)+      &
#   ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#   endif
     &                                  FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=AVERAGE(ng)%avgVwind(i,j)+      &
#   ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#   endif
     &                                  FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#   ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#   endif
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      AVERAGE(ng)%avgUwindE,                      &
     &                      AVERAGE(ng)%avgVwindN)
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
#  ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#  endif
     &                                FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)+          &
#  ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#  endif
     &                                FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=AVERAGE(ng)%avgsssflx(i,j)+    &
     &                                FORCES(ng)%sssflx(i,j,isalt)
            END DO
          END DO
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)+        &
#    ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#    endif
     &                                 FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)+        &
#    ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#    endif
     &                                 FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=AVERAGE(ng)%avguice(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%umask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=AVERAGE(ng)%avgvice(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%vmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF

        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      ICE(ng) % ui(:,:,Kout),                     &
     &                      ICE(ng) % vi(:,:,Kout),                     &
     &                      AVERAGE(ng)%avgUiceE,                       &
     &                      AVERAGE(ng)%avgViceN)
        END IF

        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=AVERAGE(ng)%avgaice(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghice(i,j)=AVERAGE(ng)%avghice(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghsno(i,j)=AVERAGE(ng)%avghsno(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtice(i,j)=AVERAGE(ng)%avgtice(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtimid(i,j)=AVERAGE(ng)%avgtimid(i,j)+      &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
#  ifdef MELT_PONDS
        IF (Aout(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgapond(i,j)=AVERAGE(ng)%avgapond(i,j)+      &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%apond(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghpond(i,j)=AVERAGE(ng)%avghpond(i,j)+      &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%hpond(i,j,Iout)
            END DO
          END DO
        END IF
#  endif
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgageice(i,j)=AVERAGE(ng)%avgageice(i,j)+    &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgiomflx(i,j)=AVERAGE(ng)%avgiomflx(i,j)+    &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig11(i,j)=AVERAGE(ng)%avgsig11(i,j)+      &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig12(i,j)=AVERAGE(ng)%avgsig12(i,j)+      &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig22(i,j)=AVERAGE(ng)%avgsig22(i,j)+      &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgT0mk(i,j)=AVERAGE(ng)%avgT0mk(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgS0mk(i,j)=AVERAGE(ng)%avgS0mk(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWfr(i,j)=AVERAGE(ng)%avgWfr(i,j)+          &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWai(i,j)=AVERAGE(ng)%avgWai(i,j)+          &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWao(i,j)=AVERAGE(ng)%avgWao(i,j)+          &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWio(i,j)=AVERAGE(ng)%avgWio(i,j)+          &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWro(i,j)=AVERAGE(ng)%avgWro(i,j)+          &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdiv,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWdiv(i,j)=AVERAGE(ng)%avgWdiv(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%wdiv(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=AVERAGE(ng)%avgutau_iw(i,j)+  &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=AVERAGE(ng)%avgchu_iw(i,j)+    &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef ICE_BIO
#  if defined BERING_10K
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgIPhL(i,j)=AVERAGE(ng)%avgIPhL(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%IcePhL(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINO3(i,j)=AVERAGE(ng)%avgINO3(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%IceNO3(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINH4(i,j)=AVERAGE(ng)%avgINH4(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%IceNH4(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgILog(i,j)=AVERAGE(ng)%avgILog(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 ICE(ng)%IceLog(i,j,Iout)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BEST_NPZ
!--------------------------------------------
! values are accumulated in the bestnpz.h file
! so do not sum them here
!--------------------------------------------
#  ifdef STATIONARY
         DO itrc=1,NTS(ng)
           DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgst(i,j,k,itrc)=                         &
#  ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#  endif
     &                     OCEAN(ng)%st(i,j,k,Nout,itrc)
               END DO
             END DO
           END DO
         END DO
# endif
# ifdef PROD3
         DO itrc=1,NPT3(ng)
           DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgpt3(i,j,k,itrc)=                        &
#  ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#  endif
     &                     OCEAN(ng)%pt3(i,j,k,Nout,itrc)
               END DO
             END DO
           END DO
         END DO
# endif
# ifdef PROD2
         DO itrc=1,NPT2(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgpt2(i,j,itrc)=                          &
#  ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#  endif
     &                     OCEAN(ng)%pt2(i,j,Nout,itrc)
               END DO
             END DO
         END DO
# endif
# ifdef STATIONARY2
         DO itrc=1,NTS2(ng)
!          DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgst2(i,j,itrc)=                          &
#  ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#  endif
     &                     OCEAN(ng)%st2(i,j,Nout,itrc)
               END DO
             END DO
!          END DO
         END DO
# endif
# ifdef BENTHIC
          DO itrc=1,NBeT(ng)
            DO k=1,NBL(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgbt(i,j,k,itrc)=                        &
     &                   AVERAGE(ng)%avgbt(i,j,k,itrc)+                 &
#  ifdef WET_DRY
     &                   GRID(ng)%rmask_full(i,j)*                      &
#  endif
     &                   OCEAN(ng)%bt(i,j,k,Nout,itrc)
                END DO
              END DO
            END DO
          END DO

# endif
# if defined CLIM_ICE_1D
          DO itrc=1,NIceT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgit(i,j,itrc)=                            &
     &                           AVERAGE(ng)%avgit(i,j,itrc)+           &
#  ifdef WET_DRY
     &                           GRID(ng)%rmask_full(i,j)*              &
#  endif
     &                           OCEAN(ng)%it(i,j,Nout,itrc)

              END DO
            END DO
          END DO
#  endif
#  ifdef BIOFLUX
         DO itrc=1,NT(ng)
           DO itrc2=1,NT(ng)
            AVERAGE(ng)%avgbflx(itrc,itrc2)=OCEAN(ng)%bflx(itrc,itrc2)
           END DO
         END DO
#  endif
# endif

# ifdef WEC
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%umask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%vmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=AVERAGE(ng)%avgu2rs(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%umask_full(i,j)*        &
#  endif
     &                                 MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=AVERAGE(ng)%avgv2rs(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%vmask_full(i,j)*        &
#  endif
     &                                 MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WEC_MELLOR
        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=AVERAGE(ng)%avgSxx2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Sxx_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=AVERAGE(ng)%avgSxy2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Sxy_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=AVERAGE(ng)%avgSyy2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Syy_bar(i,j)
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=AVERAGE(ng)%avgSxx3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Sxx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=AVERAGE(ng)%avgSxy3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Sxy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=AVERAGE(ng)%avgSyy3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Syy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=AVERAGE(ng)%avgSzx3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Szx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=AVERAGE(ng)%avgSzy3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Szy(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
# ifdef WEC
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3Sd(i,j,k)=AVERAGE(ng)%avgw3Sd(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%rmask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%W_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3St(i,j,k)=AVERAGE(ng)%avgw3St(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%rmask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%wstvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=AVERAGE(ng)%avgu3rs(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#  endif
     &                                     MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3rs(i,j,k)=AVERAGE(ng)%avgv3rs(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#  endif
     &                                     MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
# ifdef WEC_VF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=AVERAGE(ng)%avgWztw(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%zetaw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWqsp(i,j)=AVERAGE(ng)%avgWqsp(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%qsp(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbeh(i,j)=AVERAGE(ng)%avgWbeh(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%bh(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_HEIGHT
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=AVERAGE(ng)%avgWamp(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Hwave(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=AVERAGE(ng)%avgWlen(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Lwave(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (Aout(idWlep,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlep(i,j)=AVERAGE(ng)%avgWlep(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Lwavep(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_DIR
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=AVERAGE(ng)%avgWdir(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Dwave(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=AVERAGE(ng)%avgWptp(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Pwave_top(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=AVERAGE(ng)%avgWpbt(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Pwave_bot(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef BBL_MODEL
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=AVERAGE(ng)%avgWorb(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Uwave_rms(i,j)
            END DO
          END DO
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=AVERAGE(ng)%avgWdif(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Dissip_fric(i,j)
            END DO
          END DO
        END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=AVERAGE(ng)%avgWdib(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Dissip_break(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=AVERAGE(ng)%avgWdiw(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Dissip_wcap(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef ROLLER_SVENDSEN
        IF (Aout(idWbrk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbrk(i,j)=AVERAGE(ng)%avgWbrk(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Wave_break(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef WEC_ROLLER
        IF (Aout(idWdis,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdis(i,j)=AVERAGE(ng)%avgWdis(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%Dissip_roller(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWrol,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWrol(i,j)=AVERAGE(ng)%avgWrol(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 FORCES(ng)%rollA(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef UV_KIRBY
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=AVERAGE(ng)%avgUwav(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%uWave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=AVERAGE(ng)%avgVwav(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%vWave(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Accumulate vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)+    &
# ifdef WET_DRY
     &                                   GRID(ng)%pmask_full(i,j)*      &
# endif
     &                                   potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)+    &
# ifdef WET_DRY
     &                                   GRID(ng)%pmask_full(i,j)*      &
# endif
     &                                   relvor_bar(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)+  &
#  ifdef WET_DRY
     &                                       GRID(ng)%pmask_full(i,j)*  &
#  endif
     &                                       potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)+  &
#  ifdef WET_DRY
     &                                       GRID(ng)%pmask_full(i,j)*  &
#  endif
     &                                       relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Accumulate quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%rmask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%umask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%vmask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%umask_full(i,j)*      &
#  endif
     &                                   OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%vmask_full(i,j)*      &
#  endif
     &                                   OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%rmask_full(i,j)*      &
#  endif
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#  endif
     &                                     GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#  endif
     &                                     GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%rmask_full(i,j)* &
#  endif
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%umask_full(i,j)* &
#  endif
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%vmask_full(i,j)* &
#  endif
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT(i,&
     &                                                       j,k,it)+   &
#  ifdef WET_DRY
     &                                         GRID(ng)%umask_full(i,j)*&
#  endif
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT(i,&
     &                                                       j,k,it)+   &
#  ifdef WET_DRY
     &                                         GRID(ng)%vmask_full(i,j)*&
#  endif
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
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
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG(ng).eq.1) THEN
            AVGtime(ng)=time(ng)
          ELSE
            AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
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
            pfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%pmask_avg(i,j))
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_avg(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_avg(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            pfac(i,j)=fac
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
# endif
!
!  Process state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgzeta)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgzeta)
# endif
          END IF
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2d)
# endif
          END IF
        END IF

        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2d)
# endif
          END IF
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgu2dE(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2dE(i,j)
              AVERAGE(ng)%avgv2dN(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2dN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2dE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2dN)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2dE,                    &
     &                          AVERAGE(ng)%avgv2dN)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3d)
#  endif
          END IF
        END IF

        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3d)
#  endif
          END IF
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgu3dE(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3dE(i,j,k)
                AVERAGE(ng)%avgv3dN(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3dN(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3dE)
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3dN)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3dE,                    &
     &                          AVERAGE(ng)%avgv3dN)
#  endif
          END IF
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgw3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgw3d)
#  endif
          END IF
        END IF

        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgwvel(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgwvel)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgwvel)
#  endif
          END IF
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrho)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrho)
#  endif
          END IF
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=rfac(i,j)*                 &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgt(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgt(:,:,:,it))
#  endif
            END IF
          END IF
        END DO

#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=ufac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldu(:,:,it))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldu(:,:,it))
#   endif
            END IF
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=vfac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldv(:,:,it))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldv(:,:,it))
#   endif
            END IF
          END IF
        END DO
#  endif

#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=rfac(i,j)*                     &
     &               AVERAGE(ng)%avgPONsed(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgPONsed)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgPONsed)
#   endif
          END IF
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgOPALsed(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgOPALsed)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgOPALsed)
#   endif
          END IF
        END IF
        IF (Aout(idDENITsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDENITsed(i,j)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgDENITsed(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgDENITsed)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgDENITsed)
#   endif
          END IF
        END IF
        IF (Aout(idPONbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONbur(i,j)=rfac(i,j)*                     &
     &               AVERAGE(ng)%avgPONbur(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgPONbur)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgPONbur)
#   endif
          END IF
        END IF
        IF (Aout(idOPALbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALbur(i,j)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgOPALbur(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgOPALbur)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgOPALbur)
#   endif
          END IF
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
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgNPP)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgNPP)
#   endif
          END IF
        END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKv)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKv)
#   endif
          END IF
        END IF

        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKt)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKt)
#   endif
          END IF
        END IF

        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKs(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKs)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKs)
#   endif
          END IF
        END IF
#  endif

#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avghsbl(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghsbl)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghsbl)
#   endif
          END IF
        END IF
#  endif

#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avghbbl(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghbbl)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghbbl)
#   endif
          END IF
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Process 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg1)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg1)
#  endif
          END IF
        END IF

        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg2)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg2)
#  endif
          END IF
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg1)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg1)
#  endif
          END IF
        END IF

        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg2)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg2)
#  endif
          END IF
        END IF
# endif
!
!  Process surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsus)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsus)
# endif
          END IF
        END IF

        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsvs)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsvs)
# endif
          END IF
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbus)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbus)
# endif
          END IF
        END IF

        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbvs)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbvs)
# endif
          END IF
        END IF
# ifdef BBL
        IF (Aout(idUbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbrs(i,j)=rfac(i,j)*                        &
     &                                 AVERAGE(ng)%avgUbrs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUbrs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUbrs)
#  endif
          END IF
        END IF
        IF (Aout(idVbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbrs(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgVbrs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVbrs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVbrs)
#  endif
          END IF
        END IF
        IF (Aout(idUbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbws(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgUbws(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUbws)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUbws)
#  endif
          END IF
        END IF
        IF (Aout(idVbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbws(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgVbws(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVbws)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVbws)
#  endif
          END IF
        END IF
        IF (Aout(idUbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbcs(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgUbcs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUbcs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUbcs)
#  endif
          END IF
        END IF
        IF (Aout(idVbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbcs(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgVbcs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVbcs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVbcs)
#  endif
          END IF
        END IF
        IF (Aout(idUVwc,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUVwc(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgUVwc(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUVwc)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUVwc)
#  endif
          END IF
        END IF
        IF (Aout(idUbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbot(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgUbot(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUbot)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUbot)
#  endif
          END IF
        END IF
        IF (Aout(idVbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbot(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgVbot(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVbot)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVbot)
#  endif
          END IF
        END IF
        IF (Aout(idUbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUbur(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgUbur(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUbur)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUbur)
#  endif
          END IF
        END IF
        IF (Aout(idVbvr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVbvr(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgVbvr(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVbvr)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVbvr)
#  endif
          END IF
        END IF
# endif
# ifdef SOLVE3D
#  if defined RUNOFF
        IF (Aout(idRunoff,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgRunoff(i,j)=rfac(i,j)*                     &
     &                                 AVERAGE(ng)%avgRunoff(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgRunoff)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgRunoff)
#   endif
          END IF
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgPair(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgPair)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgPair)
#   endif
          END IF
        END IF
#  endif

#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgUwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwind)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwind)
#   endif
          END IF
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgVwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwind)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVwind)
#   endif
          END IF
        END IF

        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgUwindE(i,j)=rfac(i,j)*                     &
     &                                 AVERAGE(ng)%avgUwindE(i,j)
              AVERAGE(ng)%avgVwindN(i,j)=rfac(i,j)*                     &
     &                                 AVERAGE(ng)%avgVwindN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwindE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwindN)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwindE,                  &
     &                          AVERAGE(ng)%avgVwindN)
#   endif
          END IF
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgstf)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgstf)
#  endif
          END IF
        END IF

        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgswf)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgswf)
#  endif
          END IF
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=rfac(i,j)*                     &
     &                  AVERAGE(ng)%avgsssflx(i,j)
            END DO
          END DO
        END IF
#  endif

#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsrf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsrf)
#   endif
          END IF
        END IF
#  endif

#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglhf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglhf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglhf)
#   endif
          END IF
        END IF

        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgshf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgshf)
#   endif
          END IF
        END IF

        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglrf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglrf)
#   endif
          END IF
        END IF

#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgevap)
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgevap)
#    endif
          END IF
        END IF

        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgrain(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrain)
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrain)
#    endif
          END IF
        END IF
#   endif
#  endif

#  ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=ufac(i,j)*                       &
     &                AVERAGE(ng)%avguice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avguice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avguice)
#   endif
          END IF
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=vfac(i,j)*                       &
     &                    AVERAGE(ng)%avgvice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgvice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgvice)
#   endif
          END IF
        END IF

        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avguiceE(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE(ng)%avguiceE(i,j)
              AVERAGE(ng)%avgviceN(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE(ng)%avgviceN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avguiceE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgviceN)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avguiceE,                   &
     &                          AVERAGE(ng)%avgviceN)
#   endif
          END IF
        END IF

        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=rfac(i,j)*                       &
     &                       AVERAGE(ng)%avgaice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgaice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgaice)
#   endif
          END IF
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=rfac(i,j)*                       &
     &                          AVERAGE(ng)%avghice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghice)
#   endif
          END IF
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=rfac(i,j)*                       &
     &                        AVERAGE(ng)%avghsno(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghsno)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghsno)
#   endif
          END IF
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=rfac(i,j)*                       &
     &                       AVERAGE(ng)%avgtice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgtice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgtice)
#   endif
          END IF
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=rfac(i,j)*                      &
     &                    AVERAGE(ng)%avgtimid(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgtimid)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgtimid)
#   endif
          END IF
        END IF
#  ifdef MELT_PONDS
        IF (Aout(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgapond(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgapond(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgapond)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgapond)
#   endif
          END IF
        END IF
        IF (Aout(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghpond(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avghpond(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghpond)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghpond)
#   endif
          END IF
        END IF
#  endif
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=rfac(i,j)*                     &
     &                     AVERAGE(ng)%avgageice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgageice)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgageice)
#   endif
          END IF
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=rfac(i,j)*                     &
     &                     AVERAGE(ng)%avgiomflx(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgiomflx)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgiomflx)
#   endif
          END IF
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=rfac(i,j)*                      &
     &                    AVERAGE(ng)%avgsig11(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsig11)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsig11)
#   endif
          END IF
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgsig12(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsig12)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsig12)
#   endif
          END IF
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgsig22(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsig22)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsig22)
#   endif
          END IF
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgT0mk(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgT0mk)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgT0mk)
#   endif
          END IF
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgS0mk(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgS0mk)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgS0mk)
#   endif
          END IF
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWfr(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWfr)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWfr)
#   endif
          END IF
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWai(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWai)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWai)
#   endif
          END IF
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWao(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWao)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWao)
#   endif
          END IF
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWio(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWio)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWio)
#   endif
          END IF
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWro(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWro)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWro)
#  endif
          END IF
        END IF
        IF (Aout(idWdiv,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiv(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgWdiv(i,j)
             END DO
           END DO
           IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
             CALL exchange_r2d_tile (ng, tile,                          &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdiv)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdiv)
#   endif
          END IF
        END IF

        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=rfac(i,j)*                    &
     &                      AVERAGE(ng)%avgutau_iw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgutau_iw)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgutau_iw)
#  endif
          END IF
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=rfac(i,j)*                     &
     &                   AVERAGE(ng)%avgchu_iw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgchu_iw)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgchu_iw)
#   endif
          END IF
        END IF
#  endif
#  ifdef ICE_BIO
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgIPhL(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgIPhL(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgIPhL)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgIPhL)
#   endif
          END IF
        END IF
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINO3(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgINO3(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgINO3)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgINO3)
#   endif
          END IF
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINH4(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgINH4(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgINH4)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgINH4)
#   endif
          END IF
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgILog(i,j)=rfac(i,j)*AVERAGE(ng)%avgILog(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgILog)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgILog)
#   endif
          END IF
        END IF
#  endif

#  ifdef BEST_NPZ
#   ifdef BENTHIC
        DO itrc=1,NBeT(ng)
          DO k=1,NBL(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgbt(i,j,k,itrc)=                          &
     &                     rfac(i,j)*AVERAGE(ng)%avgbt(i,j,k,itrc)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, NBL(ng),     &
     &                              AVERAGE(ng)%avgbt(:,:,:,itrc))
#    ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbt(:,:,:,itrc))
#    endif
          END IF
        END DO
#   endif
#   ifdef CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgit(i,j,itrc)=                              &
     &                     rfac(i,j)*AVERAGE(ng)%avgit(i,j,itrc)
            END DO
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, NIceT(ng),     &
     &                            AVERAGE(ng)%avgit)
#    ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, NIceT(ng),         &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgit)
#    endif
          END IF
#   endif
#   ifdef STATIONARY
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                     rfac(i,j)*AVERAGE(ng)%avgst(i,j,k,itrc)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgst(:,:,:,itrc))
#    ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgst(:,:,:,itrc))
#    endif
          END IF
        END DO
#   endif
#   ifdef STATIONARY2
        DO itrc=1,NTS2(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgst2(i,j,itrc)=                             &
     &                     rfac(i,j)*AVERAGE(ng)%avgst2(i,j,itrc)
            END DO
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, NTS2(ng),      &
     &                            AVERAGE(ng)%avgst2)
#    ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, NTS2(ng),          &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgst2)
#    endif
        END IF
#   endif
#  endif
# endif

# ifdef WEC
!
!  Process radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2Sd)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2Sd)
#  endif
          END IF
        END IF

        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2Sd)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2Sd)
#  endif
          END IF
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2rs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2rs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2rs)
#  endif
          END IF
        END IF

        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2rs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2rs)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2rs)
#  endif
          END IF
        END IF
# endif
# ifdef WEC_MELLOR
        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSxx2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSxx2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxx2d)
#  endif
          END IF
        END IF

        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSxy2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSxy2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxy2d)
#  endif
          END IF
        END IF

        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSyy2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSyy2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSyy2d)
#  endif
          END IF
        END IF
# endif
# ifdef WEC
#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3Sd)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3Sd)
#   endif
          END IF
        END IF

        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3Sd)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3Sd)
#   endif
          END IF
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3rs(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3rs)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3rs)
#   endif
          END IF
        END IF

        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3RS(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3RS)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3RS)
#   endif
          END IF
        END IF

        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgW3Sd(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgW3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgW3Sd)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgW3Sd)
#   endif
          END IF
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgW3St(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgW3St(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgW3St)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgW3St)
#   endif
          END IF
        END IF
#  endif
# endif
# ifdef WEC_MELLOR
#  ifdef SOLVE3D
        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSxx3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSxx3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxx3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSxy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSxy3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxy3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSyy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSyy3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSyy3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSzx3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSzx3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSzx3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSzy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSzy3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSzy3d)
#   endif
          END IF
        END IF
#  endif
# endif
# ifdef WEC_VF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgWztw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWztw)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWztw)
#  endif
          END IF
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWqsp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWqsp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWqsp)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWqsp)
#  endif
          END IF
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbeh(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWbeh(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWbeh)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWbeh)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_HEIGHT
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWamp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWamp)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWamp)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWlen(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWlen)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWlen)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (Aout(idWlep,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlep(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWlep(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWlep)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWlep)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_DIR
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdir(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdir)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdir)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWptp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWptp)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWptp)
#  endif
          END IF
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWpbt(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWpbt)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWpbt)
#  endif
          END IF
        END IF
# endif
# ifdef BBL_MODEL
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWorb(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWorb)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWorb)
#  endif
          END IF
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdif(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdif)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdif)
#  endif
          END IF
        END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdib(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdib)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdib)
#  endif
          END IF
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdiw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdiw)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdiw)
#  endif
          END IF
        END IF
# endif
# ifdef ROLLER_SVENDSEN
        IF (Aout(idWbrk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbrk(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWbrk(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWbrk)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWbrk)
#  endif
          END IF
        END IF
# endif
# ifdef WEC_ROLLER
        IF (Aout(idWdis,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdis(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdis(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdis)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdis)
#  endif
          END IF
        END IF
        IF (Aout(idWrol,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWrol(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWrol(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWrol)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWrol)
#  endif
          END IF
        END IF
# endif
# ifdef UV_KIRBY
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgUwav(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwav)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwav)
#  endif
          END IF
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgVwav(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwav)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVwav)
#  endif
          END IF
        END IF
# endif
!
!  Process vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgpvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgpvor2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor2d)
# endif
          END IF
        END IF

        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgrvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrvor2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor2d)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgpvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgpvor3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor3d)
#  endif
          END IF
        END IF

        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgrvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrvor3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor3d)
#  endif
          END IF
        END IF
# endif
!
!  Process quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=rfac(i,j)*                         &
     &                               AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgZZ)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgZZ)
# endif
          END IF
        END IF

        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=ufac(i,j)*                         &
     &                               AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgU2)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgU2)
# endif
          END IF
        END IF

        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=vfac(i,j)*                         &
     &                               AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgV2)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgV2)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=ufac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUU)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUU)
#  endif
          END IF
        END IF

        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=vfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgVV)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVV)
#  endif
          END IF
        END IF

        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=rfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUV)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUV)
#  endif
          END IF
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHuon(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHuon)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHuon)
#  endif
          END IF
        END IF

        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHvom(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHvom)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHvom)
#  endif
          END IF
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=rfac(i,j)*                &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgTT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgTT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=ufac(i,j)*                &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgUT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgUT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=vfac(i,j)*                &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgVT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgVT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=ufac(i,j)*             &
     &                                    AVERAGE(ng)%avgHuonT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHuonT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHuonT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=vfac(i,j)*             &
     &                                    AVERAGE(ng)%avgHvomT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHvomT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHvomT(:,:,:,it))
#  endif
            END IF
          END IF
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE set_avg_tile

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!***********************************************************************
      SUBROUTINE set_detide_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            NTC, Kout,                            &
#  ifdef SOLVE3D
     &                            Nout,                                 &
#  endif
     &                            CosOmega, SinOmega,                   &
     &                            CosW_avg, CosW_sum,                   &
     &                            SinW_avg, SinW_sum,                   &
     &                            CosWCosW, SinWSinW, SinWCosW)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_ocean
      USE mod_scalars
      USE mod_tides, only: TIDES
# ifdef WET_DRY
      USE mod_grid
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
#  ifdef SOLVE3D
      integer, intent(in) :: Nout
#  endif
      integer, intent(in) :: NTC
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosOmega(:)
      real(r8), intent(in) :: SinOmega(:)
      real(r8), intent(inout) :: CosW_avg(:)
      real(r8), intent(inout) :: CosW_sum(:)
      real(r8), intent(inout) :: SinW_avg(:)
      real(r8), intent(inout) :: SinW_sum(:)
      real(r8), intent(inout) :: CosWCosW(:,:)
      real(r8), intent(inout) :: SinWSinW(:,:)
      real(r8), intent(inout) :: SinWCosW(:,:)
#  else
      real(r8), intent(in) :: CosOmega(NTC)
      real(r8), intent(in) :: SinOmega(NTC)
      real(r8), intent(inout) :: CosW_avg(NTC)
      real(r8), intent(inout) :: CosW_sum(NTC)
      real(r8), intent(inout) :: SinW_avg(NTC)
      real(r8), intent(inout) :: SinW_sum(NTC)
      real(r8), intent(inout) :: CosWCosW(NTC,NTC)
      real(r8), intent(inout) :: SinWSinW(NTC,NTC)
      real(r8), intent(inout) :: SinWCosW(NTC,NTC)
#  endif
!
!  Local variable declarations.
!
      integer :: i, it, j, k
      integer :: NTC2, mk, nk

      integer, dimension(2*NTC+1) :: indx

      real(r8) :: fac, fac1
      real(r8) :: Hsum, d

      real(r8), dimension(0:2*NTC) :: Ak
      real(r8), dimension(0:2*NTC) :: tide_harmonics
      real(r8), dimension(0:2*NTC,0:2*NTC) :: C, Y

      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

#  include "set_bounds.h"
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
      NTC2=2*NTC

      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Compute least-squares coefficients to detide time-averaged fields.
!  Notice that the coefficients are always accumulated and not
!  re-initialized.  This allows better tidal fit as the simulation
!  progresses.
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinOmega(nk)
            CosW_avg(nk)=CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Initialize.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
                TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%umask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%vmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%vbar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
#   ifdef WET_DRY
     &                                        GRID(ng)%rmask_full(i,j)* &
#   endif
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumukate Detide least-squares coefficients. They only vary in time
!  since omega (as computed in set_tides) uses model time coordinate.
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinW_avg(nk)+SinOmega(nk)
            CosW_avg(nk)=CosW_avg(nk)+CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Accumulate.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
                TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%umask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%vmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%vbar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
#   ifdef WET_DRY
     &                                         GRID(ng)%umask_full(i,j)*&
#   endif
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
!
!  Compute detide least-squares coefficients.  Build coefficient squared
!  matrix C(0:2*NTC,0:2*NTC) to invert. It is 2*NTC because we are
!  solving for real components Ak and Bk. The zero rows and column is
!  for the coefficients associated with the time mean.
!
!        F(t) = Fmean + SUM [ Ak sin(omega(k)*t) ]
!                     + SUM [ Bk cos(omega(k)*t) ]   for k=1:NTC
!
!  In the code below, all the arrays are collapsed into a single
!  dimension index such that:
!
!            k=0               mean term
!            k=1:NTC           sine terms
!            k=NTC+1:2*NTC     cosine terms
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          C(0,0)=1.0_r8                     ! time-averaged coefficient
          fac1=1.0_r8/REAL(Hcount(ng),r8)   ! global summation factor
          DO nk=1,NTC
            C(0,nk    )=fac1*SinW_sum(nk)
            C(0,nk+NTC)=fac1*CosW_sum(nk)
            C(nk,0    )=C(0,nk)             ! symmetric
            C(nk+NTC,0)=C(0,nk+NTC)         ! symmetric
            DO mk=1,NTC
              C(mk,nk)=fac1*SinWSinW(mk,nk)
              C(mk,nk+NTC)=fac1*SinWCosW(mk,nk)
              C(mk+NTC,nk)=fac1*SinWCosW(nk,mk)
              C(mk+NTC,nk+NTC)=fac1*CosWCosW(mk,nk)
            END DO
          END DO
          DO nk=0,NTC2
            DO mk=0,NTC2
              C(nk,mk)=C(mk,nk)
            END DO
          END DO
!
!  Invert least-squares coefficient matrix by LU decomposition.
!
          DO mk=0,NTC2
            DO nk=0,NTC2
              Y(mk,nk)=0.0_r8
            END DO
            Y(mk,mk)=1.0_r8                ! identity matrix
          END DO
          CALL ludcmp (C(0,0), NTC2+1, NTC2+1, indx, d)
!
!  Find inverse by columns. The matrix Y will now contain the inverse
!  of the least-squares coefficient matrix C, which will have been
!  destroyed.
!
          DO nk=0,NTC2
            CALL lubksb (C(0,0), NTC2+1, NTC2+1, indx, Y(0,nk))
          END DO
!
!  Compute time-averaged harmonics for current field average window.
!
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinW_avg(nk)
            tide_harmonics(nk+NTC)=cosW_avg(nk)
          END DO
!
!  Scale inverse by the global summation factor.
!
          DO nk=0,NTC2
            DO mk=0,NTC2
              Y(nk,mk)=fac1*Y(nk,mk)
            END DO
          END DO
        END IF
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
# ifdef WET_DRY
!  In wetting and drying, the sums are devided by the number of times
!  that each qrid point is wet. This was already computed in above in
!  routine "set_avg_tile" for nAVG steps time-average window.
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
!  Process accumulated detided averages. Notice that the regular
!  time-averaged field values used here are the ones computed in
!  routine "set_avg_tile".
!
        IF (Aout(idFsuD,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%zeta_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%zeta_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*rfac(i,j)
              END DO
              TIDES(ng)%zeta_detided(i,j)=AVERAGE(ng)%avgzeta(i,j)-     &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%ubar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%ubar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*ufac(i,j)
              END DO
              TIDES(ng)%ubar_detided(i,j)=AVERAGE(ng)%avgu2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%vbar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%vbar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*vfac(i,j)
              END DO
              TIDES(ng)%vbar_detided(i,j)=AVERAGE(ng)%avgv2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%u_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%u_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*ufac(i,j)
                END DO
                TIDES(ng)%u_detided(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%v_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%v_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*vfac(i,j)
                END DO
                TIDES(ng)%v_detided(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  Hsum=0.0_r8
                  DO nk=0,NTC2
                    Ak(nk)=0.0_r8
                    DO mk=0,NTC2
!!                    Ak(nk)=Ak(nk)+                                    &
!!   &                       Y(nk,mk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                      Ak(nk)=Ak(nk)+                                    &
     &                       Y(mk,nk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                    END DO
                    Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*rfac(i,j)
                  END DO
                  TIDES(ng)%t_detided(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k, &
     &                                                           it)-   &
     &                                          Hsum
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
      END IF

      RETURN
      END SUBROUTINE set_detide_tile
# endif

#endif
      END MODULE set_avg_mod
