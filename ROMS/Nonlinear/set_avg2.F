#include "cppdefs.h"
      MODULE set_avg2_mod
#if defined AVERAGES2
!
!svn $Id$
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
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg2

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg2 (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef WET_DRY
      USE mod_grid
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_stepping
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
      CALL set_avg2_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   KOUT,                                          &
# ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
# endif
     &                   NOUT)

# ifdef WET_DRY
      CALL set_avg_masks (ng, tile, iNLM,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % pmask_avg2,                        &
     &                    GRID(ng) % rmask_avg2,                        &
     &                    GRID(ng) % umask_avg2,                        &
     &                    GRID(ng) % vmask_avg2)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg2
!
!***********************************************************************
      SUBROUTINE set_avg2_tile (ng, tile,                               &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kout,                                    &
# ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
# endif
     &                         Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_average2
      USE mod_forces
      USE mod_grid
# ifdef SOLVE3D
      USE mod_mixing
# endif
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_ocean
      USE uv_rotate_mod, ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod, ONLY : uv_rotate3d
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
# ifdef SOLVE3D
#  ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
#  endif
      integer, intent(in) :: Nout
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j

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
      IF (nAVG2(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG2(ng)).and.                                &
     &     (MOD(iic(ng)-1,nAVG2(ng)).eq.1)).or.                         &
     &    ((iic(ng).ge.ntsAVG2(ng)).and.(nAVG2(ng).eq.1)).or.           &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN

# ifdef WET_DRY
!
!  If wetting and drying, initialize time dependent counters for wet
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg2(i,j)=MAX(0.0_r8,                        &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg2(i,j)=MAX(0.0_r8,                        &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg2(i,j)=MAX(0.0_r8,                        &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg2(i,j)=MAX(0.0_r8,                        &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
!
!  Initialize state variables.
!
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
# ifdef WET_DRY
              AVERAGE2(ng)%avgzeta(i,j)=AVERAGE2(ng)%avgzeta(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=OCEAN(ng)%ubar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE2(ng)%avgu2d(i,j)=AVERAGE2(ng)%avgu2d(i,j)*        &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=OCEAN(ng)%vbar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE2(ng)%avgv2d(i,j)=AVERAGE2(ng)%avgv2d(i,j)*        &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout2(idu2dE,ng).and.Aout2(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE2(ng)%avgu2dE,                       &
     &                      AVERAGE2(ng)%avgv2dN)
        END IF

# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu3d(i,j)=OCEAN(ng)%u(i,j,N(ng),Nout)
#  ifdef WET_DRY
              AVERAGE2(ng)%avgu3d(i,j)=AVERAGE2(ng)%avgu3d(i,j)*        &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=OCEAN(ng)%v(i,j,N(ng),Nout)
#  ifdef WET_DRY
              AVERAGE2(ng)%avgv3d(i,j)=AVERAGE2(ng)%avgv3d(i,j)*        &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout2(idu3dE,ng).and.Aout2(idv3dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % u(:,:,N(ng),Nout),              &
     &                      OCEAN(ng) % v(:,:,N(ng),Nout),              &
     &                      AVERAGE2(ng)%avgu3dE,                       &
     &                      AVERAGE2(ng)%avgv3dN)
        END IF


        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=OCEAN(ng)%rho(i,j,N(ng))
#  ifdef WET_DRY
              AVERAGE2(ng)%avgrho(i,j)=AVERAGE2(ng)%avgrho(i,j)*        &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=                            &
     &                    OCEAN(ng)%t(i,j,N(ng),Nout,itrc)
#  ifdef WET_DRY
                AVERAGE2(ng)%avgt(i,j,itrc)=AVERAGE2(ng)%avgt(i,j,it)*  &
     &                                       GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avghsbl(i,j)=AVERAGE2(ng)%avghsbl(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avghbbl(i,j)=AVERAGE2(ng)%avghbbl(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsus(i,j)=AVERAGE2(ng)%avgsus(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsvs(i,j)=AVERAGE2(ng)%avgsvs(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgbus(i,j)=AVERAGE2(ng)%avgbus(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgbvs(i,j)=AVERAGE2(ng)%avgbvs(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgstf(i,j)=AVERAGE2(ng)%avgstf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgswf(i,j)=AVERAGE2(ng)%avgswf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsrf(i,j)=AVERAGE2(ng)%avgsrf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avglhf(i,j)=AVERAGE2(ng)%avglhf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avglrf(i,j)=AVERAGE2(ng)%avglrf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgshf(i,j)=AVERAGE2(ng)%avgshf(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=FORCES(ng)%Uwind(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avguwind(i,j)=AVERAGE2(ng)%avguwind(i,j)*    &
     &                                 GRID(ng)%umask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=FORCES(ng)%Vwind(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgvwind(i,j)=AVERAGE2(ng)%avgvwind(i,j)*    &
     &                                 GRID(ng)%vmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout2(idUairE,ng).and.Aout2(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      FORCES(ng) % Uwind(:,:),                    &
     &                      FORCES(ng) % Vwind(:,:),                    &
     &                      AVERAGE2(ng)%avguwindE,                     &
     &                      AVERAGE2(ng)%avgvwindN)
        END IF

#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgevap(i,j)=AVERAGE2(ng)%avgevap(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgrain(i,j)=AVERAGE2(ng)%avgrain(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#   endif
#  endif
#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=ICE(ng)%ui(i,j,Iuout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avguice(i,j)=AVERAGE2(ng)%avguice(i,j)*      &
     &                                 GRID(ng)%umask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idUice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=ICE(ng)%vi(i,j,Iuout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgvice(i,j)=AVERAGE2(ng)%avgvice(i,j)*      &
     &                                 GRID(ng)%vmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout2(idUiceE,ng).and.Aout2(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      ICE(ng) % ui(:,:,Iuout),                    &
     &                      ICE(ng) % vi(:,:,Iuout),                    &
     &                      AVERAGE2(ng)%avguiceE,                      &
     &                      AVERAGE2(ng)%avgviceN)
        END IF


        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgaice(i,j)=ICE(ng)%ai(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgaice(i,j)=AVERAGE2(ng)%avgaice(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghice(i,j)=ICE(ng)%hi(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avghice(i,j)=AVERAGE2(ng)%avghice(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsno(i,j)=ICE(ng)%hsn(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avghsno(i,j)=AVERAGE2(ng)%avghsno(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtice(i,j)=ICE(ng)%tis(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgtice(i,j)=AVERAGE2(ng)%avgtice(i,j)*      &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtimid(i,j)=ICE(ng)%ti(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgtimid(i,j)=AVERAGE2(ng)%avgtimid(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#   ifdef MELT_PONDS
        IF (Aout2(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgapond(i,j)=ICE(ng)%apond(i,j,Iout)
#    ifdef WET_DRY
              AVERAGE2(ng)%avgapond(i,j)=AVERAGE2(ng)%avgapond(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
        IF (Aout2(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghpond(i,j)=ICE(ng)%hpond(i,j,Iout)
#    ifdef WET_DRY
              AVERAGE2(ng)%avghpond(i,j)=AVERAGE2(ng)%avghpond(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
#   endif
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgageice(i,j)=ICE(ng)%ageice(i,j,Iout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgageice(i,j)=AVERAGE2(ng)%avgageice(i,j)*  &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=ICE(ng)%io_mflux(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgiomflx(i,j)=AVERAGE2(ng)%avgiomflx(i,j)*  &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig11(i,j)=ICE(ng)%sig11(i,j,Ieout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsig11(i,j)=AVERAGE2(ng)%avgsig11(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig12(i,j)=ICE(ng)%sig12(i,j,Ieout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsig12(i,j)=AVERAGE2(ng)%avgsig12(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig22(i,j)=ICE(ng)%sig22(i,j,Ieout)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgsig22(i,j)=AVERAGE2(ng)%avgsig22(i,j)*    &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=ICE(ng)%t0mk(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgT0mk(i,j)=AVERAGE2(ng)%avgT0mk(i,j)* &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=ICE(ng)%s0mk(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgS0mk(i,j)=AVERAGE2(ng)%avgS0mk(i,j)* &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWfr(i,j)=ICE(ng)%wfr(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgWfr(i,j)=AVERAGE2(ng)%avgWfr(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWai(i,j)=ICE(ng)%wai(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgWai(i,j)=AVERAGE2(ng)%avgWai(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWao(i,j)=ICE(ng)%wao(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgWao(i,j)=AVERAGE2(ng)%avgWao(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWio(i,j)=ICE(ng)%wio(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgWio(i,j)=AVERAGE2(ng)%avgWio(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWro(i,j)=ICE(ng)%wro(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgWro(i,j)=AVERAGE2(ng)%avgWro(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=ICE(ng)%utau_iw(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgutau_iw(i,j)=AVERAGE2(ng)%avgutau_iw(i,j)*&
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=ICE(ng)%chu_iw(i,j)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgchu_iw(i,j)=AVERAGE2(ng)%avgchu_iw(i,j)*  &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef SOLVE3D
!
!  Initialized fields associated with tracers.
!
#  ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgst(i,j,itrc)=                             &
     &                       OCEAN(ng)%st(i,j,N(ng),Nout,itrc)
#   ifdef WET_DRY
              AVERAGE2(ng)%avgst(i,j,itrc)=AVERAGE2(ng)%avgst(i,j,itrc)*&
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG2(ng)) THEN
# ifdef WET_DRY
!
!  If wetting and drying, accumulate wet points counters.
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg2(i,j)=GRID(ng)%pmask_avg2(i,j)+          &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg2(i,j)=GRID(ng)%rmask_avg2(i,j)+          &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg2(i,j)=GRID(ng)%umask_avg2(i,j)+          &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg2(i,j)=GRID(ng)%vmask_avg2(i,j)+          &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
!
!  Accumulate 2D fields.
!
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=AVERAGE2(ng)%avgzeta(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=AVERAGE2(ng)%avgu2d (i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=AVERAGE2(ng)%avgv2d (i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout2(idu2dE,ng).and.Aout2(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE2(ng)%avgu2dE,                       &
     &                      AVERAGE2(ng)%avgv2dN)
        END IF


# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgu3d(i,j)=AVERAGE2(ng)%avgu3d(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%u(i,j,N(ng),Nout)
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=AVERAGE2(ng)%avgv3d(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%v(i,j,N(ng),Nout)
            END DO
          END DO
        END IF

        IF (Aout2(idu3dE,ng).and.Aout2(idv3dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % u(:,:,N(ng),Nout),              &
     &                      OCEAN(ng) % v(:,:,N(ng),Nout),              &
     &                      AVERAGE2(ng)%avgu3dE,                       &
     &                      AVERAGE2(ng)%avgv3dN)
        END IF


        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=AVERAGE2(ng)%avgrho(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%rho(i,j,N(ng))
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=AVERAGE2(ng)%avgt(i,j,itrc)+&
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  OCEAN(ng)%t(i,j,N(ng),Nout,itrc)
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=AVERAGE2(ng)%avghsbl(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=AVERAGE2(ng)%avghbbl(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=AVERAGE2(ng)%avgsus(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                 FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=AVERAGE2(ng)%avgsvs(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                 FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=AVERAGE2(ng)%avgbus(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                 FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=AVERAGE2(ng)%avgbvs(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                 FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=AVERAGE2(ng)%avgstf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=AVERAGE2(ng)%avgswf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=AVERAGE2(ng)%avgsrf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=AVERAGE2(ng)%avglhf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=AVERAGE2(ng)%avgshf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=AVERAGE2(ng)%avglrf(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=AVERAGE2(ng)%avguwind(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=AVERAGE2(ng)%avgvwind(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUairE,ng).and.Aout2(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      FORCES(ng) % Uwind(:,:),                    &
     &                      FORCES(ng) % Vwind(:,:),                    &
     &                      AVERAGE2(ng)%avguwindE,                     &
     &                      AVERAGE2(ng)%avgvwindN)
        END IF

#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=AVERAGE2(ng)%avgevap(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=AVERAGE2(ng)%avgrain(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=AVERAGE2(ng)%avguice(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%umask_full(i,j)*       &
# endif
     &                                  ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout2(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=AVERAGE2(ng)%avgvice(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%vmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF

        IF (Aout2(idUiceE,ng).and.Aout2(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      ICE(ng) % ui(:,:,Iuout),                    &
     &                      ICE(ng) % vi(:,:,Iuout),                    &
     &                      AVERAGE2(ng)%avguiceE,                      &
     &                      AVERAGE2(ng)%avgviceN)
        END IF


        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgaice(i,j)=AVERAGE2(ng)%avgaice(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghice(i,j)=AVERAGE2(ng)%avghice(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghsno(i,j)=AVERAGE2(ng)%avghsno(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtice(i,j)=AVERAGE2(ng)%avgtice(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtimid(i,j)=AVERAGE2(ng)%avgtimid(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
# ifdef MELT_PONDS
        IF (Aout2(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgapond(i,j)=AVERAGE2(ng)%avgapond(i,j)+    &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  ICE(ng)%apond(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghpond(i,j)=AVERAGE2(ng)%avghpond(i,j)+    &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  ICE(ng)%hpond(i,j,Iout)
            END DO
          END DO
        END IF
# endif
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgageice(i,j)=AVERAGE2(ng)%avgageice(i,j)+  &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=AVERAGE2(ng)%avgiomflx(i,j)+  &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig11(i,j)=AVERAGE2(ng)%avgsig11(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig12(i,j)=AVERAGE2(ng)%avgsig12(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig22(i,j)=AVERAGE2(ng)%avgsig22(i,j)+    &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=AVERAGE2(ng)%avgT0mk(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=AVERAGE2(ng)%avgS0mk(i,j)+      &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWfr(i,j)=AVERAGE2(ng)%avgWfr(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                 ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWai(i,j)=AVERAGE2(ng)%avgWai(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                 ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWao(i,j)=AVERAGE2(ng)%avgWao(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                 ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWio(i,j)=AVERAGE2(ng)%avgWio(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                 ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWro(i,j)=AVERAGE2(ng)%avgWro(i,j)+        &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                 ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=AVERAGE2(ng)%avgutau_iw(i,j)+&
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=AVERAGE2(ng)%avgchu_iw(i,j)+  &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                                  ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BIO_GOANPZ
!--------------------------------------------
! values are accumulated in the goanpz.h file
! so do not sum them here
!--------------------------------------------
        DO itrc=1,NTS(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgst(i,j,itrc)=                             &
# ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
# endif
     &                       OCEAN(ng)%st(i,j,N(ng),Nout,itrc)
            END DO
          END DO
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
      IF (((iic(ng).gt.ntsAVG2(ng)).and.                                &
     &     (MOD(iic(ng)-1,nAVG2(ng)).eq.0).and.                         &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsAVG2(ng)).and.(nAVG2(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG2(ng).eq.1) THEN
            AVG2time(ng)=time(ng)
          ELSE
            AVG2time(ng)=AVG2time(ng)+REAL(nAVG2(ng),r8)*dt(ng)
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
            pfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%pmask_avg2(i,j))
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg2(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_avg2(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_avg2(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG2(ng),r8)
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
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE2(ng)%avgzeta(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=ufac(i,j)*                      &
     &                                  AVERAGE2(ng)%avgu2d (i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=vfac(i,j)*                      &
     &                                  AVERAGE2(ng)%avgv2d (i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idu2dE,ng).and.Aout2(idv2dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE2(ng)%avgu2dE(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgu2dE(i,j)
              AVERAGE2(ng)%avgv2dN(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgv2dN(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu3d(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgu3d(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgv3d(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idu3dE,ng).and.Aout2(idv3dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE2(ng)%avgu3dE(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgu3dE(i,j)
              AVERAGE2(ng)%avgv3dN(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgv3dN(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgrho(i,j)
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=rfac(i,j)*                  &
     &                      AVERAGE2(ng)%avgt(i,j,itrc)
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avghsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avghbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Process surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgsus(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgsvs(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgbus(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgbvs(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgstf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgswf(i,j)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgsrf(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avglhf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgshf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avglrf(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE2(ng)%avguwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgvwind(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUairE,ng).and.Aout2(idVairN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE2(ng)%avguwindE(i,j)=rfac(i,j)*                    &
     &                                 AVERAGE2(ng)%avguwindE(i,j)
              AVERAGE2(ng)%avgvwindN(i,j)=rfac(i,j)*                    &
     &                                 AVERAGE2(ng)%avgvwindN(i,j)
            END DO
          END DO
        END IF

#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgevap(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgrain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif

#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE2(ng)%avguice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgvice(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUiceE,ng).and.Aout2(idViceN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE2(ng)%avguiceE(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avguiceE(i,j)
              AVERAGE2(ng)%avgviceN(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgviceN(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgaice(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgaice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghice(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avghice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsno(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avghsno(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtice(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgtice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtimid(i,j)=rfac(i,j)*                       &
     &                                 &
     &                      AVERAGE2(ng)%avgtimid(i,j)
            END DO
          END DO
        END IF
# ifdef MELT_PONDS
        IF (Aout2(idApond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgapond(i,j)=rfac(i,j)*                       &
     &                                 &
     &                      AVERAGE2(ng)%avgapond(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idHpond,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghpond(i,j)=rfac(i,j)*                       &
     &                                 &
     &                      AVERAGE2(ng)%avghpond(i,j)
            END DO
          END DO
        END IF
# endif
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgageice(i,j)=rfac(i,j)*                       &
     &                                 &
     &                      AVERAGE2(ng)%avgageice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=rfac(i,j)*                    &
     &                      AVERAGE2(ng)%avgiomflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig11(i,j)=rfac(i,j)*                     &
     &                      AVERAGE2(ng)%avgsig11(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig12(i,j)=rfac(i,j)*                     &
     &                      AVERAGE2(ng)%avgsig12(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig22(i,j)=rfac(i,j)*                     &
     &                      AVERAGE2(ng)%avgsig22(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgT0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE2(ng)%avgS0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWfr(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgWfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWai(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgWai(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWao(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgWao(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWio(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgWio(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWro(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE2(ng)%avgWro(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=rfac(i,j)*                   &
     &                      AVERAGE2(ng)%avgutau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=rfac(i,j)*                    &
     &                      AVERAGE2(ng)%avgchu_iw(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef SOLVE3D
!#  ifdef BIO_GOANPZ
!-----------------------------------------------------------
! Computes average of time series
! For production, totals are needed so this is commented out
!-----------------------------------------------------------
!        DO itrc=1,NTS(ng)
!            DO j-JstrR,JendR
!              DO i-IstrR,IendR
!                AVERAGE2(ng)%avgst(i,j,itrc)=rfac(i,j)*                &
!    &                                 AVERAGE2(ng)%avgst(i,j,itrc)
!            END DO
!          END DO
!        END DO
!#  endif
# endif
!
      END IF

      RETURN
      END SUBROUTINE set_avg2_tile
#endif
      END MODULE set_avg2_mod
