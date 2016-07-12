#include "cppdefs.h"
      MODULE ad_set_avg_mod
#if defined AD_AVERAGES && defined ADJOINT
!
!svn $Id: ad_set_avg.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  adjoint fields.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
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
      CALL wclock_on (ng, iADM, 5)
# endif
      CALL ad_set_avg_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng))

# ifdef WET_DRY
      CALL set_avg_masks (ng, tile, iADM,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % pmask_avg,                         &
     &                    GRID(ng) % rmask_avg,                         &
     &                    GRID(ng) % umask_avg,                         &
     &                    GRID(ng) % vmask_avg)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iADM, 5)
# endif
      RETURN
      END SUBROUTINE ad_set_avg
!
!***********************************************************************
      SUBROUTINE ad_set_avg_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
!
!  Local variable declarations.
!
      integer :: i, it, j, k

      real(r8) :: fac

      real(r8) :: pfac(IminS:ImaxS,JminS:JmaxS)
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
      IF (((iic(ng).lt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng),nAVG(ng)).eq.0)).or.                            &
     &    ((iic(ng).le.ntsAVG(ng)).and.(nAVG(ng).eq.1)).or.             &
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
!  Initialize adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%ad_zeta_sol(i,j)
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
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ad_ubar_sol(i,j)
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
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%ad_vbar_sol(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%ad_u(i,j,k,nstp)
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
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%ad_v(i,j,k,nstp)
#  ifdef WET_DRY
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)*    &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%ad_W(i,j,k)*        &
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

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%ad_rho(i,j,k)
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
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%ad_t(i,j,k,      &
     &                                                      nstp,it)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)*&
     &                                       GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%ad_Akv(i,j,k)
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
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%ad_Akt(i,j,k,itemp)
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
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%ad_Akt(i,j,k,isalt)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Initialize adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%ad_sustr(i,j)
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
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%ad_svstr(i,j)
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
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%ad_bustr(i,j)
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
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%ad_bvstr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%ad_stflx(i,j,itemp)
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
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%ad_stflx(i,j,isalt)
#  ifdef WET_DRY
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%ad_srflx(i,j)
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
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%ad_lhflx(i,j)
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
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%ad_lrflx(i,j)
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
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%ad_shflx(i,j)
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
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%ad_evap(i,j)
#    ifdef WET_DRY
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Initialize adjoint of quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%ad_zeta_sol(i,j)*        &
     &                               OCEAN(ng)%ad_zeta_sol(i,j)
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
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ad_ubar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_ubar_sol(i,j)
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
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%ad_vbar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_vbar_sol(i,j)
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
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%ad_u(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)
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
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%ad_v(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)
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
     &                                   (OCEAN(ng)%ad_u(i  ,j,k,nstp)+ &
     &                                    OCEAN(ng)%ad_u(i+1,j,k,nstp))*&
     &                                   (OCEAN(ng)%ad_v(i,j  ,k,nstp)+ &
     &                                    OCEAN(ng)%ad_v(i,j+1,k,nstp))
#  ifdef WET_DRY
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)*      &
     &                                   GRID(ng)%rmask_full(i,j)
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
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)*  &
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)
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
     &                                        OCEAN(ng)%ad_u(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i-1,j,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i  ,j,k,  &
     &                                                        nstp,it))
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
     &                                        OCEAN(ng)%ad_v(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i,j-1,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i,j  ,k,  &
     &                                                        nstp,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%vmask_full(i,j)
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
      ELSE IF ((iic(ng)-1).le.ntsAVG(ng)) THEN

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
!  Accumulate adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 OCEAN(ng)%ad_zeta_sol(i,j)
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
     &                                OCEAN(ng)%ad_ubar_sol(i,j)
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
     &                                OCEAN(ng)%ad_vbar_sol(i,j)
            END DO
          END DO
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
     &                                    OCEAN(ng)%ad_u(i,j,k,nstp)
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
     &                                    OCEAN(ng)%ad_v(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%ad_W(i,j,k)*        &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
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
     &                                    OCEAN(ng)%ad_rho(i,j,k)
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
     &                                       OCEAN(ng)%ad_t(i,j,k,      &
     &                                                      nstp,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%ad_Akv(i,j,k)
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
     &                                    MIXING(ng)%ad_Akt(i,j,k,itemp)
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
     &                                    MIXING(ng)%ad_Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Accumulate adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                FORCES(ng)%ad_sustr(i,j)
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
     &                                FORCES(ng)%ad_svstr(i,j)
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
     &                                FORCES(ng)%ad_bustr(i,j)
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
     &                                FORCES(ng)%ad_bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
#  ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#  endif
     &                                FORCES(ng)%ad_stflx(i,j,itemp)
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
     &                                FORCES(ng)%ad_stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%ad_srflx(i,j)
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
     &                                FORCES(ng)%ad_lhflx(i,j)
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
     &                                FORCES(ng)%ad_lrflx(i,j)
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
     &                                FORCES(ng)%ad_shflx(i,j)
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
     &                                 FORCES(ng)%ad_evap(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Accumulate adjoint quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%rmask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%ad_zeta_sol(i,j)*        &
     &                               OCEAN(ng)%ad_zeta_sol(i,j)
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
     &                               OCEAN(ng)%ad_ubar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_ubar_sol(i,j)
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
     &                               OCEAN(ng)%ad_vbar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_vbar_sol(i,j)
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
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)
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
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)
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
     &                                   (OCEAN(ng)%ad_u(i  ,j,k,nstp)+ &
     &                                    OCEAN(ng)%ad_u(i+1,j,k,nstp))*&
     &                                   (OCEAN(ng)%ad_v(i,j  ,k,nstp)+ &
     &                                    OCEAN(ng)%ad_v(i,j+1,k,nstp))
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
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)*  &
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)
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
     &                                        OCEAN(ng)%ad_u(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i-1,j,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i  ,j,k,  &
     &                                                        nstp,it))
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
     &                                        OCEAN(ng)%ad_v(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i,j-1,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i,j  ,k,  &
     &                                                        nstp,it))
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
!-----------------------------------------------------------------------
!
      IF (((iic(ng).lt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).le.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG(ng).eq.1) THEN
            AVGtime(ng)=time(ng)
          ELSE
            AVGtime(ng)=AVGtime(ng)-REAL(nAVG(ng),r8)*dt(ng)
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
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
      defined SO_SEMI
        fac=1.0_r8
#  else
        fac=1.0_r8/REAL(nAVG(ng),r8)
#  endif
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
!  Process adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
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
          END IF
        END DO

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
        END IF
#  endif
# endif
!
!  Process adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=rfac(i,j)*                        &
    &                                 AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=rfac(i,j)*                        &
    &                                 AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
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
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Process adjoint quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=rfac(i,j)*                         &
     &                               AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=ufac(i,j)*                         &
     &                               AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=vfac(i,j)*                         &
     &                               AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
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
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=rfac(i,j)*                &
     &                                        AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=ufac(i,j)*                &
     &                                        AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=vfac(i,j)*                &
     &                                        AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE ad_set_avg_tile
#endif
      END MODULE ad_set_avg_mod
