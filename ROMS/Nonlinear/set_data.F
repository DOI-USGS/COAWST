#include "cppdefs.h"
#ifdef NONLINEAR
      SUBROUTINE set_data (ng, tile)
!
!svn $Id: set_data.F 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_data_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif

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
#if defined RED_TIDE
      USE mod_biology
#endif
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
# if defined TRC_PSOURCE
      USE mod_trc_sources
# endif
!
# ifdef ANALYTICAL
      USE analytical_mod
# endif
      USE exchange_2d_mod
# ifdef DISTRIBUTE
#  if defined WET_DRY
      USE distribute_mod, ONLY : mp_boundary
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
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
      logical :: SetBC
      logical :: update = .FALSE.
# if defined WET_DRY
      logical :: bry_update
# endif

      integer :: ILB, IUB, JLB, JUB
      integer :: i, ic, itrc, j, k, my_tile

      real(r8) :: cff, cff1, cff2

# include "set_bounds.h"
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)

# ifdef NCEP_FLUXES
#  ifdef ANA_NCEP
      CALL ana_ncep_tile (ng, tile, iNLM,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng) % nustr,                        &
     &                       FORCES(ng) % nvstr,                        &
     &                       FORCES(ng) % cloud,                        &
     &                       FORCES(ng) % srflx,                        &
     &                       FORCES(ng) % lrflx,                        &
     &                       FORCES(ng) % shflx,                        &
     &                       FORCES(ng) % lhflx,                        &
     &                       FORCES(ng) % Pair,                         &
     &                       FORCES(ng) % rain,                         &
     &                       FORCES(ng) % runoff,                       &
     &                       FORCES(ng) % skt,                          &
     &                       FORCES(ng) % icec                          &
     &                   )
#  else
      CALL set_2dfld_tile (ng, tile, iNLM, idCfra,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%cloudG,                           &
     &                     FORCES(ng)%cloud,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idLrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idShea,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%shflxG,                           &
     &                     FORCES(ng)%shflx,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idLhea,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lhflxG,                           &
     &                     FORCES(ng)%lhflx,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idPair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PairG,                            &
     &                     FORCES(ng)%Pair,                             &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idRain,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%rainG,                            &
     &                     FORCES(ng)%rain,                             &
     &                     update)
#   ifdef SNOWFALL
      CALL set_2dfld_tile (ng, tile, iNLM, idsnow,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%snowG,                            &
     &                     FORCES(ng)%snow,                             &
     &                     update)
#   endif
#   ifdef SNOW_FROM_RAIN
! Subtract snowfall from total precip
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          FORCES(ng)%rain(i,j)=FORCES(ng)%rain(i,j)-                    &
     &                         FORCES(ng)%snow(i,j)
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%rain)
      END IF
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%rain)
#    endif
#   endif
#   ifdef RUNOFF
      CALL set_2dfld_tile (ng, tile, iNLM, idRunoff,                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%runoffG,                          &
     &                     FORCES(ng)%runoff,                           &
     &                     update)
#   endif
      CALL set_2dfld_tile (ng, tile, iNLM, idSkt,                       &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sktG,                             &
     &                     FORCES(ng)%skt,                              &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idIcec,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%icecG,                            &
     &                     FORCES(ng)%icec,                             &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idUnms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%nustrG,                           &
     &                     FORCES(ng)%nustr,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idVnms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%nvstrG,                           &
     &                     FORCES(ng)%nvstr,                            &
     &                     update)
#   if defined SCORRECTION || defined SRELAXATION
!
!-----------------------------------------------------------------------
!  Set surface salinity for freshwater flux correction.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, tile, iNLM, idSSSc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssG,                             &
     &                     FORCES(ng)%sss,                              &
     &                     update)
#    endif
#   endif

#   ifdef CURVGRID
!
!  If input point surface winds or interpolated from coarse data, rotate
!  to curvilinear grid.
!
      IF (.not.Linfo(1,idUsms,ng).or.                                   &
     &    (Iinfo(5,idUsms,ng).ne.Lm(ng)+2).or.                          &
     &    (Iinfo(6,idUsms,ng).ne.Mm(ng)+2)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%nustr(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%nvstr(i,j)*GRID(ng)%SinAngler(i,j)
            cff2=FORCES(ng)%nvstr(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%nustr(i,j)*GRID(ng)%SinAngler(i,j)
            FORCES(ng)%nustr(i,j)=cff1
            FORCES(ng)%nvstr(i,j)=cff2
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%nustr)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%nvstr)
        END IF
#     ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      FORCES(ng)%nustr,                           &
     &                      FORCES(ng)%nvstr)
#     endif
      END IF
#   endif
#  endif
# else

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || defined CCSM_FLUXES2D || \
    (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO_CLOUD)) && \
    !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set surface air temperature (degC).
!-----------------------------------------------------------------------
!
#   ifdef ANA_TAIR
      CALL ana_tair (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idTair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%TairG,                            &
     &                     FORCES(ng)%Tair,                             &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || defined CCSM_FLUXES2D || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO_CLOUD)) && \
     !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set surface air relative or specific humidity.
!-----------------------------------------------------------------------
!
#   ifdef ANA_HUMIDITY
      CALL ana_humid (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idQair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HairG,                            &
     &                     FORCES(ng)%Hair,                             &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef SOLVE3D

#   if defined ALBEDO_FILE && defined SHORTWAVE
      CALL set_2dfld_tile (ng, tile, iNLM, idAlbe,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%albedoG,                          &
     &                     FORCES(ng)%albedo,                           &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#    ifdef ICE_MODEL
      FORCES(ng)%albedo_ice = FORCES(ng)%albedo
#    endif
#   endif
#  if defined CLOUDS && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set cloud fraction (nondimensional).  Notice that clouds are
!  processed first in case that they are used to adjust shortwave
!  radiation.
!-----------------------------------------------------------------------
!
#    ifdef ANA_CLOUD
      CALL ana_cloud (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, tile, iNLM, idCfra,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%cloudG,                           &
     &                     FORCES(ng)%cloud,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif

#  if defined SHORTWAVE && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SRFLUX
      CALL ana_srflux (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef DIURNAL_SRFLUX
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      CALL ana_srflux (ng, tile, iNLM)
#   endif
#  endif
#  if defined RED_TIDE && defined DAILY_SHORTWAVE
!
!-----------------------------------------------------------------------
!  Set kinematic daily-averaged surface solar shortwave radiation flux
!  (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idAsrf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG_avg,                       &
     &                     FORCES(ng)%srflx_avg,                        &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ANA_LRFLUX
      CALL ana_lrflux (ng, tile, iNLM)
#  elif defined BULK_FLUXES && !defined LONGWAVE && !defined LONGWAVE_OUT \
   && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idLrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined LONGWAVE_OUT && defined BULK_FLUXES && (!defined ANA_LRFLUX \
   && !defined AIR_OCEAN)
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idLdwn,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

#  if (defined BULK_FLUXES || defined BULK_FLUXES2D || defined ECOSIM \
       || defined SPECTRAL_LIGHT) && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set surface winds (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_WINDS
      CALL ana_winds (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idUair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%UwindG,                           &
     &                     FORCES(ng)%Uwind,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfld_tile (ng, tile, iNLM, idVair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%VwindG,                           &
     &                     FORCES(ng)%Vwind,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN

#    ifdef CURVGRID
!
!  If input point surface winds or interpolated from coarse data, rotate
!  to curvilinear grid.
!
      IF (.not.Linfo(1,idUair,ng).or.                                   &
     &    (Iinfo(5,idUair,ng).ne.Lm(ng)+2).or.                          &
     &    (Iinfo(6,idUair,ng).ne.Mm(ng)+2)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%Uwind(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%Vwind(i,j)*GRID(ng)%SinAngler(i,j)
            cff2=FORCES(ng)%Vwind(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%Uwind(i,j)*GRID(ng)%SinAngler(i,j)
            FORCES(ng)%Uwind(i,j)=cff1
            FORCES(ng)%Vwind(i,j)=cff2
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%UWind)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%VWind)
        END IF

#     ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      FORCES(ng)%UWind,                           &
     &                      FORCES(ng)%VWind)
#     endif
      END IF
#    endif
#   endif
#  endif

# ifdef SOLVE3D
#  if defined BULK_FLUXES && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set rain fall rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_RAIN
      CALL ana_rain (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idrain,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%rainG,                            &
     &                     FORCES(ng)%rain,                             &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#   ifdef SNOWFALL
#    ifdef ANA_SNOW
      CALL ana_snow (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, tile, iNLM, idsnow,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%snowG,                            &
     &                     FORCES(ng)%snow,                             &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
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
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if !defined BULK_FLUXES && !defined ATM2OCN_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_STFLUX
      CALL ana_stflux (ng, tile, iNLM, itemp)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%stflx (:,:,itemp),                &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
       defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, set surface net heat flux from first NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itemp),             &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stflxG(:,:,:,itemp),            &
     &                       FORCES(ng)%stflx (:,:,itemp),              &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Set sea surface temperature (SST) and heat flux sensitivity to
!  SST (dQdSST) which are used for surface heat flux correction.
!-----------------------------------------------------------------------
!
#   ifdef ANA_SST
      CALL ana_sst (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idSSTc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sstG,                             &
     &                     FORCES(ng)%sst,                              &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
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
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_BTFLUX
      CALL ana_btflux (ng, tile, iNLM, itemp)
#  else
      CALL set_2dfld_tile (ng, tile, iNLM, idTbot(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%btflx (:,:,itemp),                &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SALINITY
#   if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
        defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, set kinematic surface freshwater (E-P) flux
!  (m/s) from first NLM run.
!-----------------------------------------------------------------------
!
#   else
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater (E-P) flux (m/s).
!-----------------------------------------------------------------------
!
#   endif
#   ifdef ANA_SSFLUX
      CALL ana_stflux (ng, tile, iNLM, isalt)
#   else
#    if defined NL_BULK_FLUXES && !defined BULK_FLUXES
      CALL set_2dfld_tile (ng, tile, iNLM, idEmPf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%stflx (:,:,isalt),                &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN

#    elif !(defined EMINUSP || defined SRELAXATION)
      CALL set_2dfld_tile (ng, tile, iNLM, idsfwf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%stflx (:,:,isalt),                &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif

#   if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
        defined BULK_FLUXES && defined NL_BULK_FLUXES
      IF (Nrun.gt.1) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idEmPf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stflxG(:,:,:,isalt),            &
     &                       FORCES(ng)%stflx (:,:,isalt),              &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#   endif

#   if defined SCORRECTION || defined SRELAXATION
!
!-----------------------------------------------------------------------
!  Set surface salinity for freshwater flux correction.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, tile, iNLM, idSSSc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssG,                             &
     &                     FORCES(ng)%sss,                              &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_BSFLUX
      CALL ana_btflux (ng, tile, iNLM, isalt)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idTbot(isalt),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%btflx (:,:,isalt),                &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
!
!-----------------------------------------------------------------------
!  Set kinematic surface and bottom pasive tracer fluxes (T m/s).
!-----------------------------------------------------------------------
!
      DO itrc=NAT+1,NT(ng)
#   ifdef ANA_SPFLUX
        CALL ana_stflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%stflx (:,:,itrc),               &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#   ifdef ANA_BPFLUX
        CALL ana_btflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idTbot(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%btflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%btflx (:,:,itrc),               &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#   endif
      END DO
#  endif
# endif

# ifndef AIR_OCEAN
#  if !defined BULK_FLUXES && !defined BULK_FLUXES2D
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SMFLUX
      CALL ana_smflux (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, tile, iNLM, idUsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sustrG,                           &
     &                     FORCES(ng)%sustr,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfld_tile (ng, tile, iNLM, idVsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%svstrG,                           &
     &                     FORCES(ng)%svstr,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN

#    ifdef CURVGRID
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

#     ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      FORCES(ng)%sustr,                           &
     &                      FORCES(ng)%svstr)
#     endif
      END IF
#    endif
#   endif
#  endif
# endif
# endif

# if !defined AIR_OCEAN   && defined FOUR_DVAR      && \
      defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, set surface wind stress components from first
!  NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idUsms,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%sustrG,                         &
     &                       FORCES(ng)%sustr,                          &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_2dfld_tile (ng, tile, iNLM, idVsms,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%svstrG,                         &
     &                       FORCES(ng)%svstr,                          &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

#  if (defined BULK_FLUXES || defined ECOSIM || \
       defined SPECTRAL_LIGHT || defined CCSM_FLUXES2D ||  \
       defined ATM_PRESS) && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Set surface air pressure (mb).
!-----------------------------------------------------------------------
!
#  ifdef ANA_PAIR
      CALL ana_pair (ng, tile, iNLM)
#  else
      SetBC=.TRUE.
!     SetBC=.FALSE.
      CALL set_2dfld_tile (ng, tile, iNLM, idPair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PairG,                            &
     &                     FORCES(ng)%Pair,                             &
     &                     update, SetBC)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# ifdef CICE_MODEL
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!  Here, R/Cp = 0.286
          FORCES(ng)%PotT(i,j)=FORCES(ng)%Tair(i,j)*                    &
     &          (1000.0/FORCES(ng)%Pair(i,j))**0.286
        END DO
      END DO
# endif
# if defined WAVE_DATA && !defined INWAVE_MODEL
!
!-----------------------------------------------------------------------
!  Set surface wind-induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
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
#   endif

#   ifdef WAVES_HEIGHT
      CALL set_2dfld_tile (ng, tile, iNLM, idWamp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HwaveG,                           &
     &                     FORCES(ng)%Hwave,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef WAVES_LENGTH
      CALL set_2dfld_tile (ng, tile, iNLM, idWlen,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%LwaveG,                           &
     &                     FORCES(ng)%Lwave,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#   ifdef WAVES_LENGTHP
      CALL set_2dfld_tile (ng, tile, iNLM, idWlep,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%LwavepG,                          &
     &                     FORCES(ng)%Lwavep,                           &
     &                     update)
#   endif
#   ifdef WAVES_TOP_PERIOD
      CALL set_2dfld_tile (ng, tile, iNLM, idWptp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_topG,                       &
     &                     FORCES(ng)%Pwave_top,                        &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef WAVES_BOT_PERIOD
      CALL set_2dfld_tile (ng, tile, iNLM, idWpbt,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_botG,                       &
     &                     FORCES(ng)%Pwave_bot,                        &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined WAVES_UB
      CALL set_2dfld_tile (ng, tile, iNLM, idWorb,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Uwave_rmsG,                       &
     &                     FORCES(ng)%Uwave_rms,                        &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#   if defined WAVES_DISS
      CALL set_2dfld_tile (ng, tile, iNLM, idWdib,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Dissip_breakG,                    &
     &                     FORCES(ng)%Dissip_break,                     &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iNLM, idWdiw,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Dissip_wcapG,                     &
     &                     FORCES(ng)%Dissip_wcap,                      &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#   if defined ROLLER_SVENDSEN
      CALL set_2dfld_tile (ng, tile, iNLM, idWbrk,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Wave_breakG,                      &
     &                     FORCES(ng)%Wave_break,                       &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
# endif

# if defined SSSFLX && defined SOLVE3D
      CALL set_2dfld_tile (ng, tile, iNLM, idSSSf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssflxG,                          &
     &                     FORCES(ng)%sssflx,                           &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined FASTICE_CLIMATOLOGY && defined SOLVE3D
      CALL set_2dfld_tile (ng, tile, iNLM, idFastIce,                   &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%fastice_clmG,                     &
     &                     FORCES(ng)%fastice_clm,                      &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined SOLVE3D && (defined ECOSIM || defined SPECTRAL_LIGHT)
!
!-----------------------------------------------------------------------
!  Compute spectral irradiance and cosine of average zenith angle of
!  downwelling spectral photons.
!-----------------------------------------------------------------------
!
      CALL ana_specir (ng, tile, iNLM)
# endif

# ifdef ANA_SPINNING
!
!-----------------------------------------------------------------------
!  Set time-varying rotation force (centripetal accelerations) for
!  polar coordinate grids.
!-----------------------------------------------------------------------
!
      CALL ana_spinning (ng, tile, iNLM)
# endif

# if defined TRC_PSOURCE
#  ifdef ANA_TRC_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks for passive tracers.
!-----------------------------------------------------------------------
!
      CALL ana_trc_psource (ng, tile, iNLM)
#  endif
# endif

!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
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
          IF (exit_flag.ne.NoError) RETURN

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
        DO itrc=1,NT(ng)
          IF (LtracerSrc(itrc,ng)) THEN
#   ifdef ONE_TRACER_SOURCE
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, 1, 1,            &
     &                      1, 1, 1,                                    &
     &                      SOURCES(ng) % TsrcG(:,itrc),                &
     &                      SOURCES(ng) % Tsrc(itrc),                   &
     &                      update)
#   elif defined TWO_D_TRACER_SOURCE
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Nsrc(ng), 1,     &
     &                      1, Nsrc(ng), 1,                             &
     &                      SOURCES(ng) % TsrcG(:,:,itrc),              &
     &                      SOURCES(ng) % Tsrc(:,itrc),                 &
     &                      update)
#   else
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Nsrc(ng), N(ng), &
     &                      1, Nsrc(ng), N(ng),                         &
     &                      SOURCES(ng) % TsrcG(:,:,:,itrc),            &
     &                      SOURCES(ng) % Tsrc(:,:,itrc),               &
     &                      update)
#   endif
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
      END IF
# endif

# ifdef RED_TIDE
!
!-----------------------------------------------------------------------
!  Red tide Observed Dissolved Inorganic Nutrient.
!-----------------------------------------------------------------------
!
      CALL set_3dfld_tile (ng, tile, iNLM, idODIN,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng)%DIN_obsG,                          &
     &                     OCEAN(ng)%DIN_obs,                           &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.  In grid refinement, only the
!  coarser grid (RefineScale(ng)=0) open boundary conditions data is
!  processed and needed.
!-----------------------------------------------------------------------
!
!  Free-surface
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
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
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(ieast,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_east,                  &
     &                      BOUNDARY(ng) % zeta_east,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(isouth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % zetaG_south,                 &
     &                      BOUNDARY(ng) % zeta_south,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(inorth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_north,                 &
     &                      BOUNDARY(ng) % zeta_north,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif

# if defined WET_DRY
!
! Ensure that water level on boundary cells is above bed elevation.
!
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrR,JendR
              cff=Dcrit(ng)-GRID(ng)%h(0,j)
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

        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrR,JendR
              cff=Dcrit(ng)-GRID(ng)%h(Lm(ng)+1,j)
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

        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrR,IendR
              cff=Dcrit(ng)-GRID(ng)%h(i,0)
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

        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrR,IendR
              cff=Dcrit(ng)-GRID(ng)%h(i,Mm(ng)+1)
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
!  2D momentum.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
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
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(iwest,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(iwest), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_west,                  &
     &                      BOUNDARY(ng) % vbar_west,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(ieast,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_east,                  &
     &                      BOUNDARY(ng) % ubar_east,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(ieast,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(ieast), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_east,                  &
     &                      BOUNDARY(ng) % vbar_east,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(isouth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(isouth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_south,                 &
     &                      BOUNDARY(ng) % ubar_south,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(isouth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_south,                 &
     &                      BOUNDARY(ng) % vbar_south,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(inorth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(inorth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_north,                 &
     &                      BOUNDARY(ng) % ubar_north,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(inorth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_north,                 &
     &                      BOUNDARY(ng) % vbar_north,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

#ifdef ICE_MODEL
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_AIOBC
        CALL ana_aiobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isAice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idAibc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % aiG_west,                    &
     &                      BOUNDARY(ng) % ai_west,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isAice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idAibc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % aiG_east,                    &
     &                      BOUNDARY(ng) % ai_east,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isAice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idAibc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % aiG_south,                   &
     &                      BOUNDARY(ng) % ai_south,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isAice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idAibc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % aiG_north,                   &
     &                      BOUNDARY(ng) % ai_north,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_HIOBC
        CALL ana_hiobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isHice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHibc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hiG_west,                    &
     &                      BOUNDARY(ng) % hi_west,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isHice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHibc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hiG_east,                    &
     &                      BOUNDARY(ng) % hi_east,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isHice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHibc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % hiG_south,                   &
     &                      BOUNDARY(ng) % hi_south,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isHice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHibc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hiG_north,                   &
     &                      BOUNDARY(ng) % hi_north,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# if defined ICE_BIO
#  ifdef BERING_10K
        CALL ana_IcePhLbc (ng, tile, iNLM)
        CALL ana_IceNO3bc (ng, tile, iNLM)
        CALL ana_IceNH4bc (ng, tile, iNLM)
#  endif
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_HSNOBC
        CALL ana_hsnobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isHsno,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHsnbc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hsnG_west,                   &
     &                      BOUNDARY(ng) % hsn_west,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isHsno,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHsnbc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hsnG_east,                   &
     &                      BOUNDARY(ng) % hsn_east,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isHsno,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHsnbc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % hsnG_south,                  &
     &                      BOUNDARY(ng) % hsn_south,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isHsno,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHsnbc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hsnG_north,                  &
     &                      BOUNDARY(ng) % hsn_north,                   &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_TIOBC
        CALL ana_tiobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isTice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idTibc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % tiG_west,                    &
     &                      BOUNDARY(ng) % ti_west,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isTice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idTibc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % tiG_east,                    &
     &                      BOUNDARY(ng) % ti_east,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isTice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idTibc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % tiG_south,                   &
     &                      BOUNDARY(ng) % ti_south,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isTice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idTibc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % tiG_north,                   &
     &                      BOUNDARY(ng) % ti_north,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef MELT_PONDS
#  ifdef ANA_SFWATOBC
        CALL ana_sfwatobc (ng, tile, iNLM)
#  else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isApond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idApdbc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % apondG_west,                 &
     &                      BOUNDARY(ng) % apond_west,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isApond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idApdbc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % apondG_east,                 &
     &                      BOUNDARY(ng) % apond_east,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isApond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idApdbc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % apondG_south,                &
     &                      BOUNDARY(ng) % apond_south,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isApond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idApdbc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % apondG_north,                &
     &                      BOUNDARY(ng) % apond_north,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(iwest,isHpond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHpdbc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hpondG_west,                 &
     &                      BOUNDARY(ng) % hpond_west,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isHpond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHpdbc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hpondG_east,                 &
     &                      BOUNDARY(ng) % hpond_east,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isHpond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHpdbc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % hpondG_south,                &
     &                      BOUNDARY(ng) % hpond_south,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isHpond,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idHpdbc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % hpondG_north,                &
     &                      BOUNDARY(ng) % hpond_north,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
#  endif
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_SIG11OBC
        CALL ana_sig11obc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isSig11,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS11bc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig11G_west,                 &
     &                      BOUNDARY(ng) % sig11_west,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isSig11,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS11bc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig11G_east,                 &
     &                      BOUNDARY(ng) % sig11_east,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isSig11,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS11bc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % sig11G_south,                &
     &                      BOUNDARY(ng) % sig11_south,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isSig11,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS11bc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig11G_north,                &
     &                      BOUNDARY(ng) % sig11_north,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_SIG22OBC
        CALL ana_sig22obc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isSig22,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS22bc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig22G_west,                 &
     &                      BOUNDARY(ng) % sig22_west,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isSig22,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS22bc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig22G_east,                 &
     &                      BOUNDARY(ng) % sig22_east,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isSig22,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS22bc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % sig22G_south,                &
     &                      BOUNDARY(ng) % sig22_south,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isSig22,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS22bc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig22G_north,                &
     &                      BOUNDARY(ng) % sig22_north,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_SIG12OBC
        CALL ana_sig12obc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isSig12,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS12bc(iwest), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig12G_west,                 &
     &                      BOUNDARY(ng) % sig12_west,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isSig12,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS12bc(ieast), JLB, JUB, 1,      &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig12G_east,                 &
     &                      BOUNDARY(ng) % sig12_east,                  &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isSig12,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS12bc(isouth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % sig12G_south,                &
     &                      BOUNDARY(ng) % sig12_south,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isSig12,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idS12bc(inorth), ILB, IUB, 1,     &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % sig12G_north,                &
     &                      BOUNDARY(ng) % sig12_north,                 &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_MIOBC
        CALL ana_miobc (ng, tile, iNLM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isUice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idUibc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % uiG_west,                    &
     &                      BOUNDARY(ng) % ui_west,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(iwest,isVice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idVibc(iwest), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % viG_west,                    &
     &                      BOUNDARY(ng) % vi_west,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isUice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idUibc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % uiG_east,                    &
     &                      BOUNDARY(ng) % ui_east,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(ieast,isVice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idVibc(ieast), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % viG_east,                    &
     &                      BOUNDARY(ng) % vi_east,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isUice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idUibc(isouth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % uiG_south,                   &
     &                      BOUNDARY(ng) % ui_south,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(isouth,isVice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idVibc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % viG_south,                   &
     &                      BOUNDARY(ng) % vi_south,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isUice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idUibc(inorth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % uiG_north,                   &
     &                      BOUNDARY(ng) % ui_north,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF (LBC(inorth,isVice,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idVibc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % viG_north,                   &
     &                      BOUNDARY(ng) % vi_north,                    &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF
#endif

# ifdef SOLVE3D
!
!  3D momentum.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
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
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(iwest,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(iwest), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_west,                     &
     &                      BOUNDARY(ng) % v_west,                      &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(ieast,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(ieast), JLB, JUB, N(ng),   &
     &                      0, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_east,                     &
     &                      BOUNDARY(ng) % u_east,                      &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(ieast,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(ieast), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_east,                     &
     &                      BOUNDARY(ng) % v_east,                      &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(isouth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(isouth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_south,                    &
     &                      BOUNDARY(ng) % u_south,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(isouth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(isouth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_south,                    &
     &                      BOUNDARY(ng) % v_south,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(inorth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(inorth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_north,                    &
     &                      BOUNDARY(ng) % u_north,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (LBC(inorth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(inorth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_north,                    &
     &                      BOUNDARY(ng) % v_north,                     &
     &                      update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
#  endif
      END IF
!
!  Tracers.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
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
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (LBC(ieast,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(ieast,itrc),             &
     &                        JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_east(:,:,:,itrc),       &
     &                        BOUNDARY(ng) % t_east(:,:,itrc),          &
     &                        update)
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (LBC(isouth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(isouth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_south(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_south(:,:,itrc),         &
     &                        update)
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (LBC(inorth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(inorth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_north(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_north(:,:,itrc),         &
     &                        update)
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
        END IF
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set sea surface height climatology (m).
!-----------------------------------------------------------------------
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
        IF (exit_flag.ne.NoError) RETURN
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
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
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
        IF (exit_flag.ne.NoError) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%vbarclmG,                        &
     &                       CLIMA(ng)%vbarclm,                         &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
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
!-----------------------------------------------------------------------
!  Set 3D momentum climatology (m/s).
!-----------------------------------------------------------------------
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
        IF (exit_flag.ne.NoError) RETURN
!
        CALL set_3dfld_tile (ng, tile, iNLM, idVclm,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%vclmG,                           &
     &                       CLIMA(ng)%vclm,                            &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
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
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng)%u(:,:,:,i),OCEAN(ng)%v(:,:,:,i))
      END DO
#   endif

#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Set tracers climatology.
!-----------------------------------------------------------------------
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
          IF (exit_flag.ne.NoError) RETURN
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

#  ifdef AICLIMATOLOGY
#   ifdef ANA_AICLIMA
        CALL ana_aiclima (ng, tile, iNLM)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idAIclm,                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%aiclmG,                          &
     &                       CLIMA(ng)%aiclm,                           &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
#  ifdef HICLIMATOLOGY
#   ifdef ANA_HICLIMA
        CALL ana_hiclima (ng, tile, iNLM)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idHIclm,                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%hiclmG,                          &
     &                       CLIMA(ng)%hiclm,                           &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef MICLIMATOLOGY
#   ifdef ANA_MICLIMA
        CALL ana_miclima (ng, tile, iNLM)
#   else
        CALL set_2dfld_tile (ng, tile, iNLM, idUIclm,                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%uiclmG,                          &
     &                       CLIMA(ng)%uiclm,                           &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
        CALL set_2dfld_tile (ng, tile, iNLM, idVIclm,                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%viclmG,                          &
     &                       CLIMA(ng)%viclm,                           &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
# endif

# if defined SOLVE3D && defined OCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set omega climatology.
!-----------------------------------------------------------------------
!
#  ifdef ANA_OCLIMA
      CALL ana_oclima (ng, tile, iNLM)
#  else
      CALL set_3dfld_tile (ng, tile, iNLM, idOclm,                      &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       CLIMA(ng)%oclmG,                           &
     &                       CLIMA(ng)%oclm,                            &
     &                       update)
#  endif
#  if defined OFFLINE
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            OCEAN(ng)%W(i,j,k) = CLIMA(ng)%oclm(i,j,k)/                 &
     &          ( (GRID(ng)%pm(i,j)) * (GRID(ng)%pn(i,j)) )
          ENDDO
        ENDDO
      ENDDO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          OCEAN(ng)%W(:,:,:))
      END IF
#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng)%W)
#   endif
#  endif
# endif

# if defined SOLVE3D && defined AKTCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set vertical diffusion climatology.
!-----------------------------------------------------------------------
!
#  ifdef ANA_AKTCLIMA
      CALL ana_aktclima (ng, tile, iNLM)
#  else
      CALL set_3dfld_tile (ng, tile, iNLM, idAclm,                      &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       CLIMA(ng)%AclmG,                           &
     &                       CLIMA(ng)%Aclm,                            &
     &                       update)
#  endif
#  if defined OFFLINE
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
#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng)%Akt(:,:,:,2))
#   endif
#  endif
# endif

# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Set weak contraint forcing.
!-----------------------------------------------------------------------
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
        IF (exit_flag.ne.NoError) RETURN

#  ifndef SOLVE3D
!
!  Set 2D momentum forcing.
!
        CALL set_2dfld_tile (ng, tile, iNLM, idUbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%ubarG,                           &
     &                       OCEAN(ng)%f_ubar,                          &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_2dfld_tile (ng, tile, iNLM, idVbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%vbarG,                           &
     &                       OCEAN(ng)%f_vbar,                          &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SOLVE3D
!
!  Set 3D momentum.
!
        CALL set_3dfld_tile (ng, tile, iNLM, idUvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%uG,                              &
     &                       OCEAN(ng)%f_u,                             &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_3dfld_tile (ng, tile, iNLM, idVvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%vG,                              &
     &                       OCEAN(ng)%f_v,                             &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
!
!  Set 3D tracers.
!
        DO itrc=1,NT(ng)
          CALL set_3dfld_tile (ng, tile, iNLM, idTvar(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OCEAN(ng)%tG(:,:,:,:,itrc),              &
     &                         OCEAN(ng)%f_t(:,:,:,itrc),               &
     &                         update)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#  endif
      END IF
# endif

      RETURN
      END SUBROUTINE set_data_tile
#else
      SUBROUTINE set_data
      RETURN
      END SUBROUTINE set_data
#endif
