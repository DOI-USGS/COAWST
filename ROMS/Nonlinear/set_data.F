#include "cppdefs.h"
#ifdef NONLINEAR
      SUBROUTINE set_data (ng, tile)
!
!svn $Id: set_data.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  other input data. It time-interpolates between snapshots.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4, __LINE__, MyFile)
# endif
      CALL set_data_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE set_data
!
!***********************************************************************
      SUBROUTINE set_data_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
# if defined HYPOXIA_SRM || defined RED_TIDE
      USE mod_biology
# endif
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
      USE mod_sources
# if defined INWAVE_MODEL && defined INWAVE_SWAN_COUPLING
      USE bndwavez_mod
      USE bndwaveuv_mod
# endif

# if defined TRC_PSOURCE
      USE mod_trc_sources
# endif
!
# ifdef ANALYTICAL
      USE analytical_mod
# endif
      USE exchange_2d_mod
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
# ifdef DISTRIBUTE
#  if defined WET_DRY
      USE distribute_mod,  ONLY : mp_boundary
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE strings_mod,     ONLY : FoundError
# if defined OFFLINE
      USE mod_stepping
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
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
!
!  Local variable declarations.
!
      logical :: Lprocess, SetBC
      logical :: update = .FALSE.
# if defined WET_DRY
      logical :: bry_update
# endif
!
      integer :: ILB, IUB, JLB, JUB
      integer :: i, ic, itrc, j, k, my_tile
!
      real(r8) :: cff, cff1, cff2
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", set_data_tile"
!
# include "set_bounds.h"
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!=======================================================================
!  Set point Sources/Sinks (river runoff).
!=======================================================================
!
!  Point Source/Sink vertically integrated mass transport.
!
# ifdef ANA_PSOURCE
      IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
        CALL ana_psource (ng, tile, iNLM)
      END IF
# else
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (LuvSrc(ng).or.LwSrc(ng)) THEN
          CALL set_ngfld (ng, iNLM, idRtra, 1, Nsrc(ng), 1,             &
     &                    1, Nsrc(ng), 1,                               &
     &                    SOURCES(ng) % QbarG,                          &
     &                    SOURCES(ng) % Qbar,                           &
     &                    update)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=1,Nsrc(ng)
              SOURCES(ng)%Qsrc(i,k)=SOURCES(ng)%Qbar(i)*                &
     &                              SOURCES(ng)%Qshape(i,k)
            END DO
          END DO
#  endif
        END IF

#  ifdef SOLVE3D
!
!  Tracer Sources/Sinks.
!
        DO itrc=1,NT(ng)
          IF (LtracerSrc(itrc,ng)) THEN
#    ifdef ONE_TRACER_SOURCE
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, 1, 1,            &
     &                      1, 1, 1,                                    &
     &                      SOURCES(ng) % TsrcG(:,itrc),                &
     &                      SOURCES(ng) % Tsrc(itrc),                   &
     &                      update)
#    elif defined TWO_D_TRACER_SOURCE
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Nsrc(ng), 1,     &
     &                      1, Nsrc(ng), 1,                             &
     &                      SOURCES(ng) % TsrcG(:,:,itrc),              &
     &                      SOURCES(ng) % Tsrc(:,itrc),                 &
     &                      update)
#    else
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Nsrc(ng), N(ng), &
     &                      1, Nsrc(ng), N(ng),                         &
     &                      SOURCES(ng) % TsrcG(:,:,:,itrc),            &
     &                      SOURCES(ng) % Tsrc(:,:,itrc),               &
     &                      update)
#    endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
      END IF
# endif
!
!=======================================================================
!  Set forcing data.
!=======================================================================
!
!  Set switch to process surface atmospheric fields.
!
# if defined FOUR_DVAR   && \
     defined BULK_FLUXES && defined PRIOR_BULK_FLUXES
!  In 4D-Var data assimilation applications, the user have the option
!  to fix the prior (background phase) surface fluxes in the successive
!  outer loops (Nouter>1) or the final analysis phase. In such case, the
!  fluxes are read from the background trajectory.
!
      IF (Nrun.eq.1) THEN
        Lprocess=.TRUE.
      ELSE
        Lprocess=.FALSE.
      END IF
# else
      Lprocess=.TRUE.
# endif

# ifdef SOLVE3D

#  ifdef CLOUDS
!
!  Set cloud fraction (nondimensional).  Notice that clouds are
!  processed first in case that they are used to adjust shortwave
!  radiation.
!
      IF (Lprocess) THEN
#   ifdef ANA_CLOUD
        CALL ana_cloud (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
        CALL set_2dfld_tile (ng, tile, iNLM, idCfra,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%cloudG,                         &
     &                       FORCES(ng)%cloud,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
      END IF
#  endif

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || \
    (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO))
!
!  Set surface air temperature (degC).
!
      IF (Lprocess) THEN
#   ifdef ANA_TAIR
        CALL ana_tair (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
        CALL set_2dfld_tile (ng, tile, iNLM, idTair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%TairG,                          &
     &                       FORCES(ng)%Tair,                           &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
      END IF
#  endif

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO))
!
!  Set surface air relative or specific humidity.
!
      IF (Lprocess) THEN
#   ifdef ANA_HUMIDITY
        CALL ana_humid (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
        CALL set_2dfld_tile (ng, tile, iNLM, idQair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%HairG,                          &
     &                       FORCES(ng)%Hair,                           &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
      END IF
#  endif

#  ifdef SHORTWAVE
!
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!
#   ifdef ANA_SRFLUX
      CALL ana_srflux (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
      CALL set_2dfld_tile (ng, tile, iNLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef DIURNAL_SRFLUX
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      IF (Lprocess) THEN
        CALL ana_srflux (ng, tile, iNLM)
      END IF
#   endif
#  endif

#  if defined RED_TIDE && defined DAILY_SHORTWAVE
!
!  Set kinematic daily-averaged surface solar shortwave radiation flux
!  (degC m/s).
!
      CALL set_2dfld_tile (ng, tile, iNLM, idAsrf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG_avg,                       &
     &                     FORCES(ng)%srflx_avg,                        &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined BULK_FLUXES     && \
    !(defined LONGWAVE        || defined LONGWAVE_OUT) && \
     ( (defined FRC_COUPLING  && defined TIME_INTERP)  || \
       !defined FRC_COUPLING )
!
!  Surface net longwave radiation (degC m/s).
!
      IF (Lprocess) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idLrad,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%lrflxG,                         &
     &                       FORCES(ng)%lrflx,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  if defined LONGWAVE_OUT   && defined BULK_FLUXES  && \
     ( (defined FRC_COUPLING && defined TIME_INTERP) || \
       !defined FRC_COUPLING )
!
!  Surface downwelling longwave radiation (degC m/s).
!
      IF (Lprocess) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idLdwn,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%lrflxG,                         &
     &                       FORCES(ng)%lrflx,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  if (defined BULK_FLUXES || defined BULK_FLUXES2D || defined ECOSIM \
       || defined SPECTRAL_LIGHT)
!
!  Set surface winds (m/s).
!
      IF (Lprocess) THEN
#   ifdef ANA_WINDS
        CALL ana_winds (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
        CALL set_2dfld_tile (ng, tile, iNLM, idUair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%UwindG,                         &
     &                       FORCES(ng)%Uwind,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%VwindG,                         &
     &                       FORCES(ng)%Vwind,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef CURVGRID
!
!  If input point surface winds or interpolated from coarse data, rotate
!  to curvilinear grid.
!
        IF (.not.Linfo(1,idUair,ng).or.                                 &
     &      (Iinfo(5,idUair,ng).ne.Lm(ng)+2).or.                        &
     &      (Iinfo(6,idUair,ng).ne.Mm(ng)+2)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff1=FORCES(ng)%Uwind(i,j)*GRID(ng)%CosAngler(i,j)+       &
     &             FORCES(ng)%Vwind(i,j)*GRID(ng)%SinAngler(i,j)
              cff2=FORCES(ng)%Vwind(i,j)*GRID(ng)%CosAngler(i,j)-       &
     &             FORCES(ng)%Uwind(i,j)*GRID(ng)%SinAngler(i,j)
              FORCES(ng)%Uwind(i,j)=cff1
              FORCES(ng)%Vwind(i,j)=cff2
            END DO
          END DO

          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              FORCES(ng)%UWind)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              FORCES(ng)%VWind)
          END IF

#     ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng)%UWind,                         &
     &                        FORCES(ng)%VWind)
#     endif
        END IF
#    endif
#   endif
      END IF
#  endif

#  ifdef BULK_FLUXES
!
!  Set rain fall rate (kg/m2/s).
!
      IF (Lprocess) THEN
#   ifdef ANA_RAIN
        CALL ana_rain (ng, tile, iNLM)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
        CALL set_2dfld_tile (ng, tile, iNLM, idrain,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%rainG,                          &
     &                       FORCES(ng)%rain,                           &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef SNOWFALL
#    ifdef ANA_SNOW
      CALL ana_snow (ng, tile, iNLM)
#    elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
      CALL set_2dfld_tile (ng, tile, iNLM, idsnow,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%snowG,                            &
     &                     FORCES(ng)%snow,                             &
     &                     update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif
!
!-----------------------------------------------------------------------
! Fresh water runoff (kg/m2/s)
!-----------------------------------------------------------------------
!
#   ifdef RUNOFF
      CALL set_2dfld_tile (ng, tile, iNLM, idRunoff,                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%runoffG,                          &
     &                     FORCES(ng)%runoff,                           &
     &                     update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
      END IF
#  endif

#  ifndef BULK_FLUXES
!
!  Set kinematic surface net heat flux (degC m/s).
!
#   ifdef ANA_STFLUX
      CALL ana_stflux (ng, tile, iNLM, itemp)
#   elif (defined FRC_COUPLING && defined TIME_INTERP) || \
         !defined FRC_COUPLING
      CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stfluxG(:,:,:,itemp),             &
     &                     FORCES(ng)%stflux (:,:,itemp),               &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef QCORRECTION
!
!  Set sea surface temperature (SST) and heat flux sensitivity to
!  SST (dQdSST) which are used for surface heat flux correction.
!
#   ifdef ANA_SST
      CALL ana_sst (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idSSTc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sstG,                             &
     &                     FORCES(ng)%sst,                              &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
!
#   ifdef ANA_DQDSST
      CALL ana_dqdsst (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, iddQdT,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%dqdtG,                            &
     &                     FORCES(ng)%dqdt,                             &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
!
!  Set kinematic bottom net heat flux (degC m/s).
!
#  ifdef ANA_BTFLUX
      CALL ana_btflux (ng, tile, iNLM, itemp)
#  else
      CALL set_2dfld_tile (ng, tile, iNLM, idTbot(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btfluxG(:,:,:,itemp),             &
     &                     FORCES(ng)%btflux (:,:,itemp),               &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SALINITY
#   ifdef ANA_SSFLUX
!
!  Surface freshwater (E-P) flux (m/s) from analytical function.
!
      CALL ana_stflux (ng, tile, iNLM, isalt)
#   else

#    if !(defined BULK_FLUXES  || defined EMINUSP      || \
          defined FRC_COUPLING)
!
!  Surface freshwater (E-P) flux (m/s) from NetCDF variable "swflux".
!
      CALL set_2dfld_tile (ng, tile, iNLM, idsfwf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stfluxG(:,:,:,isalt),             &
     &                     FORCES(ng)%stflux (:,:,isalt),               &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    elif defined BULK_FLUXES  && !defined EMINUSP
!
!  Surface freshwater (E-P) flux (m/s) from NetCDF variable "EminusP".
!
      IF (Lprocess) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idEmPf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stfluxG(:,:,:,isalt),           &
     &                       FORCES(ng)%stflux (:,:,isalt),             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#    endif
#   endif

#   if defined SCORRECTION || defined SRELAXATION
!
!  Set surface salinity for freshwater flux correction.
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, tile, iNLM, idSSSc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssG,                             &
     &                     FORCES(ng)%sss,                              &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif
!
!  Set kinematic bottom salt flux (m/s).
!
#   ifdef ANA_BSFLUX
      CALL ana_btflux (ng, tile, iNLM, isalt)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idTbot(isalt),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btfluxG(:,:,:,isalt),             &
     &                     FORCES(ng)%btflux (:,:,isalt),               &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
!
!  Set kinematic surface and bottom passive tracer fluxes (T m/s).
!
      DO itrc=NAT+1,NT(ng)
#   ifdef ANA_SPFLUX
        CALL ana_stflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stfluxG(:,:,:,itrc),            &
     &                       FORCES(ng)%stflux (:,:,itrc),              &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef ANA_BPFLUX
        CALL ana_btflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idTbot(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%btfluxG(:,:,:,itrc),            &
     &                       FORCES(ng)%btflux (:,:,itrc),              &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
      END DO
#  endif
# endif

#  if !defined BULK_FLUXES && !defined BULK_FLUXES2D
!
!  Set kinematic surface momentum flux (m2/s2).
!
#  ifdef ANA_SMFLUX
      CALL ana_smflux (ng, tile, iNLM)
#  elif (defined FRC_COUPLING && defined TIME_INTERP) || \
        !defined FRC_COUPLING
      CALL set_2dfld_tile (ng, tile, iNLM, idUsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sustrG,                           &
     &                     FORCES(ng)%sustr,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL set_2dfld_tile (ng, tile, iNLM, idVsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%svstrG,                           &
     &                     FORCES(ng)%svstr,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef CURVGRID
!
!  If input point wind stress, rotate to curvilinear grid.  Notice
!  that rotation is done at RHO-points. It does not matter.
!
      IF (.not.Linfo(1,idUsms,ng).or.                                   &
     &    (Iinfo(5,idUsms,ng).ne.Lm(ng)+1).or.                          &
     &    (Iinfo(6,idUsms,ng).ne.Mm(ng)+2)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%sustr(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%svstr(i,j)*GRID(ng)%SinAngler(i,j)
            cff2=FORCES(ng)%svstr(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%sustr(i,j)*GRID(ng)%SinAngler(i,j)
            FORCES(ng)%sustr(i,j)=cff1
            FORCES(ng)%svstr(i,j)=cff2
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%sustr)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%svstr)
        END IF

#    ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      FORCES(ng)%sustr,                           &
     &                      FORCES(ng)%svstr)
#    endif
      END IF
#   endif
#  endif
# endif

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || defined ATM_PRESS)
!
!  Set surface air pressure (mb).
!
      IF (Lprocess) THEN
#  ifdef ANA_PAIR
        CALL ana_pair (ng, tile, iNLM)
#  elif (defined FRC_COUPLING && defined TIME_INTERP) || \
        !defined FRC_COUPLING
        SetBC=.TRUE.
!       SetBC=.FALSE.
        CALL set_2dfld_tile (ng, tile, iNLM, idPair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%PairG,                          &
     &                       FORCES(ng)%Pair,                           &
     &                       update, SetBC)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END IF
# endif
# if defined WAVE_DATA && !defined INWAVE_MODEL
!
!  Set surface wind-induced wave amplitude, direction and period.
!
#  ifdef ANA_WWAVE
      CALL ana_wwave (ng, tile, iNLM)
#  else
#   ifdef WAVES_DIR
      CALL set_2dfld_tile (ng, tile, iNLM, idWdir,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%DwaveG,                           &
     &                     FORCES(ng)%Dwave,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef WAVES_DIRP
      CALL set_2dfld_tile (ng, tile, iNLM, idWdip,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%DwavepG,                           &
     &                     FORCES(ng)%Dwavep,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef WAVES_HEIGHT
      CALL set_2dfld_tile (ng, tile, iNLM, idWamp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HwaveG,                           &
     &                     FORCES(ng)%Hwave,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif

#   ifdef WAVES_LENGTH
      CALL set_2dfld_tile (ng, tile, iNLM, idWlen,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%LwaveG,                           &
     &                     FORCES(ng)%Lwave,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif
#   ifdef WAVES_LENGTHP
      CALL set_2dfld_tile (ng, tile, iNLM, idWlep,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%LwavepG,                          &
     &                     FORCES(ng)%Lwavep,                           &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef WAVES_TOP_PERIOD
      CALL set_2dfld_tile (ng, tile, iNLM, idWptp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_topG,                       &
     &                     FORCES(ng)%Pwave_top,                        &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif

#   ifdef WAVES_BOT_PERIOD
      CALL set_2dfld_tile (ng, tile, iNLM, idWpbt,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_botG,                       &
     &                     FORCES(ng)%Pwave_bot,                        &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif

#   if defined WAVES_UB
      CALL set_2dfld_tile (ng, tile, iNLM, idWorb,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Uwave_rmsG,                       &
     &                     FORCES(ng)%Uwave_rms,                        &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif
#   if defined WAVES_DISS
      CALL set_2dfld_tile (ng, tile, iNLM, idWdib,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Dissip_breakG,                    &
     &                     FORCES(ng)%Dissip_break,                     &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL set_2dfld_tile (ng, tile, iNLM, idWdiw,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Dissip_wcapG,                     &
     &                     FORCES(ng)%Dissip_wcap,                      &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#   endif
#   if defined ROLLER_SVENDSEN
      CALL set_2dfld_tile (ng, tile, iNLM, idWbrk,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Wave_breakG,                      &
     &                     FORCES(ng)%Wave_break,                       &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
# endif

# if defined SSSFLX && defined SOLVE3D
      CALL set_2dfld_tile (ng, tile, iNLM, idSSSf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssflxG,                          &
     &                     FORCES(ng)%sssflx,                           &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined FASTICE_CLIMATOLOGY && defined SOLVE3D
      CALL set_2dfld_tile (ng, tile, iNLM, idFastIce,                   &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%fastice_clmG,                     &
     &                     FORCES(ng)%fastice_clm,                      &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SOLVE3D && (defined ECOSIM || defined SPECTRAL_LIGHT)
!
!  Compute spectral irradiance and cosine of average zenith angle of
!  downwelling spectral photons.
!
      CALL ana_specir (ng, tile, iNLM)
# endif

# ifdef ANA_SPINNING
!
!  Set time-varying rotation force (centripetal accelerations) for
!  polar coordinate grids.
!
      CALL ana_spinning (ng, tile, iNLM)
# endif

# ifdef HYPOXIA_SRM
!
!  Total respiration rate for hypoxia.
!
#  ifdef ANA_RESPIRATION
      CALL ana_respiration (ng, tile, iNLM)
#  else
      CALL set_3dfld_tile (ng, tile, iNLM, idResR,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng)%respirationG,                      &
     &                     OCEAN(ng)%respiration,                       &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef RED_TIDE
!
!  Red tide Observed Dissolved Inorganic Nutrient.
!
      CALL set_3dfld_tile (ng, tile, iNLM, idODIN,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng)%DIN_obsG,                          &
     &                     OCEAN(ng)%DIN_obs,                           &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

#if defined FOUR_DVAR   && \
    defined BULK_FLUXES && defined PRIOR_BULK_FLUXES
!
!=======================================================================
!  In 4D-Var data assimilation algorithms, the user has the option to
!  impose the prior (background phase) surface fluxes in the successive
!  outer loops (Nouter>1) or the final analysis phase. Such fluxes were
!  computed by routine "bulk_fluxes" and stored in the NLM background
!  trajectory, BLK structure. Notice that "bulk_fluxes" is called in
!  "main3d" only during the prior trajectory computation.
!
!  Therefore, the fluxes are time interpolated from the pior snapshots.
!=======================================================================
!
      IF (.not.Lprocess) THEN
!
!  Set prior surface wind stress components.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idUsms,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%sustrG,                         &
     &                       FORCES(ng)%sustr,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVsms,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%svstrG,                         &
     &                       FORCES(ng)%svstr,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ATM_PRESS
!
!  Set prior surface air pressure.
!
        SetBC=.TRUE.
!!      SetBC=.FALSE.
        CALL set_2dfld_tile (ng, tile, iNLM, idPair,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%PairG,                          &
     &                       FORCES(ng)%Pair,                           &
     &                       update, SetBC)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  ifndef ANA_STFLUX
!
!  Set prior surface net heat flux.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itemp),             &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stfluxG(:,:,:,itemp),           &
     &                       FORCES(ng)%stflux (:,:,itemp),             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SALINITY && !defined ANA_SSFLUX
!
!  Set prior surface freshwater flux.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idEmPf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stfluxG(:,:,:,isalt),           &
     &                       FORCES(ng)%stflux (:,:,isalt),             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
      END IF
#endif
!
!=======================================================================
!  Set open boundary conditions fields.
!=======================================================================
!
!  Set free-surface open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
# ifdef ANA_FSOBC
        CALL ana_fsobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_west,                  &
     &                      BOUNDARY(ng) % zeta_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(ieast,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_east,                  &
     &                      BOUNDARY(ng) % zeta_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(isouth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % zetaG_south,                 &
     &                      BOUNDARY(ng) % zeta_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(inorth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_north,                 &
     &                      BOUNDARY(ng) % zeta_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
# endif
# if defined INWAVE_MODEL && defined INWAVE_SWAN_COUPLING
        CALL bndwavez (ng, tile)
# endif
# if defined WET_DRY
!
!  Ensure that water level on boundary cells is above bed elevation.
!
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrR,JendR
!jcw          cff=Dcrit(ng)-GRID(ng)%h(0,j)
              cff=0.95_r8*Dcrit(ng)-GRID(ng)%h(0,j)
              IF (BOUNDARY(ng)%zeta_west(j).le.cff) THEN
                BOUNDARY(ng)%zeta_west(j)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_west)
#  endif
        END IF
!
        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrR,JendR
!jcw          cff=Dcrit(ng)-GRID(ng)%h(Lm(ng)+1,j)
              cff=0.95_r8*Dcrit(ng)-GRID(ng)%h(Lm(ng)+1,j)
              IF (BOUNDARY(ng)%zeta_east(j).le.cff) THEN
                BOUNDARY(ng)%zeta_east(j)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_east)
#  endif
        END IF
!
        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrR,IendR
!jcw          cff=Dcrit(ng)-GRID(ng)%h(i,0)
              cff=0.95_r8*Dcrit(ng)-GRID(ng)%h(i,0)
              IF (BOUNDARY(ng)%zeta_south(i).le.cff) THEN
                BOUNDARY(ng)%zeta_south(i)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_south)
#  endif
        END IF
!
        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrR,IendR
!jcw          cff=Dcrit(ng)-GRID(ng)%h(i,Mm(ng)+1)
              cff=0.95_r8*Dcrit(ng)-GRID(ng)%h(i,Mm(ng)+1)
              IF (BOUNDARY(ng)%zeta_north(i).le.cff) THEN
                BOUNDARY(ng)%zeta_north(i)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_north)
#  endif
        END IF
# endif
      END IF
!
!  Set 2D momentum component open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
# ifdef ANA_M2OBC
        CALL ana_m2obc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_west,                  &
     &                      BOUNDARY(ng) % ubar_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(iwest,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(iwest), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_west,                  &
     &                      BOUNDARY(ng) % vbar_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(ieast,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_east,                  &
     &                      BOUNDARY(ng) % ubar_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(ieast,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(ieast), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_east,                  &
     &                      BOUNDARY(ng) % vbar_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(isouth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(isouth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_south,                 &
     &                      BOUNDARY(ng) % ubar_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(isouth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_south,                 &
     &                      BOUNDARY(ng) % vbar_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(inorth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(inorth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_north,                 &
     &                      BOUNDARY(ng) % ubar_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(inorth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_north,                 &
     &                      BOUNDARY(ng) % vbar_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
# endif
# if defined INWAVE_MODEL && defined INWAVE_SWAN_COUPLING
        CALL bndwaveuv (ng, tile)
# endif
      END IF

# ifdef SOLVE3D
!
!  Set 3D momentum components open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
#  ifdef ANA_M3OBC
        CALL ana_m3obc (ng, tile, iNLM)
#  else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(iwest), JLB, JUB, N(ng),   &
     &                      0, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_west,                     &
     &                      BOUNDARY(ng) % u_west,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(iwest,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(iwest), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_west,                     &
     &                      BOUNDARY(ng) % v_west,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(ieast,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(ieast), JLB, JUB, N(ng),   &
     &                      0, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_east,                     &
     &                      BOUNDARY(ng) % u_east,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(ieast,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(ieast), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_east,                     &
     &                      BOUNDARY(ng) % v_east,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(isouth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(isouth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_south,                    &
     &                      BOUNDARY(ng) % u_south,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(isouth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(isouth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_south,                    &
     &                      BOUNDARY(ng) % v_south,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(inorth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(inorth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_north,                    &
     &                      BOUNDARY(ng) % u_north,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (LBC(inorth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(inorth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_north,                    &
     &                      BOUNDARY(ng) % v_north,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#  endif
      END IF
!
!  Set tracer variables open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
#  ifdef ANA_TOBC
        CALL ana_tobc (ng, tile, iNLM)
#  else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
           DO itrc=1,NT(ng)
            IF (LBC(iwest,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(iwest,itrc),             &
     &                        JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_west(:,:,:,itrc),       &
     &                        BOUNDARY(ng) % t_west(:,:,itrc),          &
     &                        update)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
!
            IF (LBC(ieast,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(ieast,itrc),             &
     &                        JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_east(:,:,:,itrc),       &
     &                        BOUNDARY(ng) % t_east(:,:,itrc),          &
     &                        update)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
!
            IF (LBC(isouth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(isouth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_south(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_south(:,:,itrc),         &
     &                        update)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
!
            IF (LBC(inorth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(inorth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_north(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_north(:,:,itrc),         &
     &                        update)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
        END IF
#  endif
      END IF
# endif
!
!=======================================================================
!  Set climatology data.
!=====================================================================
!
!  Set sea surface height climatology (m).
!
      IF (LsshCLM(ng)) THEN
# ifdef ANA_SSH
        CALL ana_ssh (ng, tile, iNLM)
# else
        CALL set_2dfld_tile (ng, tile, iNLM, idSSHc,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%sshG,                            &
     &                       CLIMA(ng)%ssh,                             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined OFFLINE
        DO j=Jstr,Jend
          DO i=Istr,Iend
            OCEAN(ng)%zeta(i,j,1) = CLIMA(ng)%ssh(i,j)
            OCEAN(ng)%zeta(i,j,2) = CLIMA(ng)%ssh(i,j)
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_r2d_tile (ng, tile,                           &
     &                      LBi, UBi, LBj, UBj, OCEAN(ng)%zeta(:,:,i))
          END DO
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%zeta)
#  endif
# endif
      END IF
!
!  Set 2D momentum components climatology (m/s).
!
      IF (Lm2CLM(ng)) THEN
# ifdef ANA_M2CLIMA
        CALL ana_m2clima (ng, tile, iNLM)
# else
        CALL set_2dfld_tile (ng, tile, iNLM, idUbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%ubarclmG,                        &
     &                       CLIMA(ng)%ubarclm,                         &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%vbarclmG,                        &
     &                       CLIMA(ng)%vbarclm,                         &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined OFFLINE
        DO j=Jstr,Jend
          DO i=Istr,Iend
            OCEAN(ng)%ubar(i,j,1) = CLIMA(ng)%ubarclm(i,j)
            OCEAN(ng)%ubar(i,j,2) = CLIMA(ng)%ubarclm(i,j)
            OCEAN(ng)%vbar(i,j,1) = CLIMA(ng)%vbarclm(i,j)
            OCEAN(ng)%vbar(i,j,2) = CLIMA(ng)%vbarclm(i,j)
          ENDDO
        ENDDO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_u2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj, OCEAN(ng)%ubar(:,:,i))
            CALL exchange_v2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj, OCEAN(ng)%vbar(:,:,i))
          END DO
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%ubar, OCEAN(ng)%vbar)
#  endif
# endif
      END IF

# ifdef SOLVE3D
!
!  Set 3D momentum components climatology (m/s).
!
      IF (Lm3CLM(ng)) THEN
#  ifdef ANA_M3CLIMA
        CALL ana_m3clima (ng, tile, iNLM)
#  else
        CALL set_3dfld_tile (ng, tile, iNLM, idUclm,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%uclmG,                           &
     &                       CLIMA(ng)%uclm,                            &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL set_3dfld_tile (ng, tile, iNLM, idVclm,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%vclmG,                           &
     &                       CLIMA(ng)%vclm,                            &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  if defined OFFLINE
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              OCEAN(ng)%u(i,j,k,nnew(ng)) = CLIMA(ng)%uclm(i,j,k)
              OCEAN(ng)%v(i,j,k,nnew(ng)) = CLIMA(ng)%vclm(i,j,k)
            END DO
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_u3d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OCEAN(ng)%u(:,:,:,i))
            CALL exchange_v3d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OCEAN(ng)%v(:,:,:,i))
          END DO
        END IF
#   ifdef DISTRIBUTE
        DO i=1,2
          CALL mp_exchange3d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng)%u(:,:,:,i),OCEAN(ng)%v(:,:,:,i))
      END DO
#   endif
#  endif
      END IF
!
!  Set tracer variables climatology.
!
#  ifdef ANA_TCLIMA
      IF (ANY(LtracerCLM(:,ng))) THEN
        CALL ana_tclima (ng, tile, iNLM)
      END IF
#  else
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng)) THEN
          ic=ic+1
          CALL set_3dfld_tile (ng, tile, iNLM, idTclm(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         CLIMA(ng)%tclmG(:,:,:,:,ic),             &
     &                         CLIMA(ng)%tclm (:,:,:,ic),               &
     &                         update)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
#  if defined OFFLINE
      DO itrc=1,NAT
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              OCEAN(ng)%t(i,j,k,nnew(ng),itrc) =                        &
     &                  CLIMA(ng)%tclm(i,j,k,itrc)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NAT
          DO i=1,2
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              OCEAN(ng)%t(:,:,:,i,itrc))
          END DO
        END DO
      END IF
#   ifdef DISTRIBUTE
      DO itrc=1,NAT
        DO i=1,2
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng)%t(:,:,:,i,itrc))
        END DO
      END DO
#   endif
#  endif

#  if defined SOLVE3D && defined AKTCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set vertical diffusion climatology.
!-----------------------------------------------------------------------
!
#   ifdef ANA_AKTCLIMA
      CALL ana_aktclima (ng, tile, iNLM)
#   else
      CALL set_3dfld_tile (ng, tile, iNLM, idAclm,                      &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       CLIMA(ng)%AclmG,                           &
     &                       CLIMA(ng)%Aclm,                            &
     &                       update)
#   endif
#   if defined OFFLINE
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            MIXING(ng)%Akt(i,j,k,2) = CLIMA(ng)%Aclm(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          MIXING(ng)%Akt(:,:,:,2))
      END IF
#    ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng)%Akt(:,:,:,2))
#    endif
#   endif
#  endif
# endif

# if defined NLM_OUTER                || \
     defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined SP4DVAR
!
!=======================================================================
!  Set weak contraint forcing.
!=======================================================================
!
      IF (FrequentImpulse(ng)) THEN
!
!  Set free-surface forcing.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idFsur,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%zetaG,                           &
     &                       OCEAN(ng)%f_zeta,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifndef SOLVE3D
!
!  Set 2D momentum forcing.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idUbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%ubarG,                           &
     &                       OCEAN(ng)%f_ubar,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL set_2dfld_tile (ng, tile, iNLM, idVbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%vbarG,                           &
     &                       OCEAN(ng)%f_vbar,                          &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  else
!
!  Set 3D momentum.
!
        CALL set_3dfld_tile (ng, tile, iNLM, idUvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%uG,                              &
     &                       OCEAN(ng)%f_u,                             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL set_3dfld_tile (ng, tile, iNLM, idVvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%vG,                              &
     &                       OCEAN(ng)%f_v,                             &
     &                       update)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set 3D tracers.
!
        DO itrc=1,NT(ng)
          CALL set_3dfld_tile (ng, tile, iNLM, idTvar(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OCEAN(ng)%tG(:,:,:,:,itrc),              &
     &                         OCEAN(ng)%f_t(:,:,:,itrc),               &
     &                         update)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
#  endif
      END IF
# endif
!
      RETURN
      END SUBROUTINE set_data_tile
#else
      SUBROUTINE set_data
      RETURN
      END SUBROUTINE set_data
#endif
