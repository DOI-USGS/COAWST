#include "cppdefs.h"
      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
#if defined RED_TIDE
      USE mod_biology
#endif
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined NLM_OUTER || defined RED_TIDE            || \
    defined TLM_CHECK || defined TL_W4DPSAS          || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
      USE mod_ocean
#endif
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
#endif
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================

#ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG(:,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

# ifdef SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), update(1),                        &
#  ifdef ONE_TRACER_SOURCE
     &                    1, 1, 1, 2, 1, 1, 1,                          &
     &                    SOURCES(ng) % TsrcG(:,i))
#  elif defined TWO_D_TRACER_SOURCE
     &                    1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,            &
     &                    SOURCES(ng) % TsrcG(:,:,i))
#  else
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
#  endif
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif
#endif

#if !defined ANA_NCEP && defined NCEP_FLUXES
!
!  Surface wind stress components.
!
      CALL get_2dfld (ng, iNLM, idUnms, ncFRCid(idUnms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % nustrG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng, iNLM, idVnms, ncFRCid(idVnms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % nvstrG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Cloud fraction.
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % cloudG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Surface solar shortwave radiation.
!
      CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % srflxG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Surface net longwave radiation.
!
      CALL get_2dfld (ng, iNLM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % lrflxG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Net sensible heat flux
!
      CALL get_2dfld (ng, iNLM, idShea, ncFRCid(idShea,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % shflxG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Net latent heat flux
!
      CALL get_2dfld (ng, iNLM, idLhea, ncFRCid(idLhea,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % lhflxG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Surface air pressure.
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % PairG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Rainfall
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % rainG)
      IF (exit_flag.ne.NoError) RETURN
# if defined SNOWFALL && !defined ANA_SNOW
!
!  Snowfall
!
      CALL get_2dfld (ng, iNLM, idsnow, ncFRCid(idsnow,ng),             &
     &                nFfiles(ng), FRCname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % snowG)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Skin temperature
!
      CALL get_2dfld (ng, iNLM, idSkt, ncFRCid(idSkt,ng),               &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % sktG)
      IF (exit_flag.ne.NoError) RETURN
!
!  Observed (NCEP) ice concentration
!
      CALL get_2dfld (ng, iNLM, idIcec, ncFRCid(idIcec,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % icecG)
      IF (exit_flag.ne.NoError) RETURN
!
# if !defined ANA_RUNOFF && defined RUNOFF
!
!-----------------------------------------------------------------------
!  Surface runoff rate
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idRunoff, ncFRCid(idRunoff,ng),         &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % runoffG)
      IF (exit_flag.ne.NoError) RETURN
# endif

#endif

#if (!defined ANA_WINDS && !defined AIR_OCEAN && \
     !defined NCEP_FLUXES) && \
    (defined BULK_FLUXES2D || defined BULK_FLUXES || \
     defined ECOSIM || defined SPECTRAL_LIGHT)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % UwindG)
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % VwindG)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if !defined AIR_OCEAN && !defined NCEP_FLUXES
# if !defined ANA_SMFLUX && !defined BULK_FLUXES && !defined ANA_NCEP \
          && !defined BULK_FLUXES2D
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FORCES(ng) % sustrG)
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iNLM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                FORCES(ng) % svstrG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined FOUR_DVAR && defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface wind stress components from initial
!  NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idUsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  FORCES(ng) % sustrG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  FORCES(ng) % svstrG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif
#endif

#ifdef SOLVE3D

# if !defined ANA_SSSFLX && defined SSSFLX
      CALL get_2dfld (ng, iNLM, idSSSf, ncFRCid(idSSSf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % sssflxG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_FASTICE && defined FASTICE_CLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idFastIce, ncFRCid(idFastIce,ng),       &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % fastice_clmG)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if (!defined ANA_PAIR && !defined AIR_OCEAN) && (defined BULK_FLUXES \
      || defined ECOSIM || defined CCSM_FLUXES2D || defined ATM_PRESS \
      || defined SPECTRAL_LIGHT)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % PairG)
      IF (exit_flag.ne.NoError) RETURN
#endif
#if !defined ANA_WWAVE && defined WAVE_DATA && !defined INWAVE_MODEL
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
# ifdef WAVES_DIR
      CALL get_2dfld (ng, iNLM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % DwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_HEIGHT
      CALL get_2dfld (ng, iNLM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % HwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_LENGTH
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % LwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef WAVES_LENGTHP
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlep,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % LwavepG)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef WAVES_TOP_PERIOD
      CALL get_2dfld (ng, iNLM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_topG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_BOT_PERIOD
      CALL get_2dfld (ng, iNLM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_botG(:,:,1))
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined WAVES_UB
      CALL get_2dfld (ng, iNLM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Uwave_rmsG)
      IF (exit_flag.ne.NoError) RETURN
# endif
# if defined WAVES_DISS
      CALL get_2dfld (ng, iNLM, idWdib, ncFRCid(idWdib,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Dissip_breakG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng, iNLM, idWdiw, ncFRCid(idWdiw,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Dissip_wcapG)
      IF (exit_flag.ne.NoError) RETURN
# endif
# if defined ROLLER_SVENDSEN
      CALL get_2dfld (ng, iNLM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Wave_breakG)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

# if (!defined ANA_TAIR && !defined AIR_OCEAN) && \
    ( defined BULK_FLUXES || defined ECOSIM || defined CCSM_FLUXES2D || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO_CLOUD) \
     || defined SPECTRAL_LIGHT)
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % TairG)
      IF (exit_flag.ne.NoError) RETURN
#endif

# if (!defined ANA_HUMIDITY && !defined AIR_OCEAN) && \
    ( defined BULK_FLUXES || defined ECOSIM || defined CCSM_FLUXES2D || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO_CLOUD) \
    || defined SPECTRAL_LIGHT)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % HairG)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef SOLVE3D

# if defined ALBEDO_FILE && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Albedo.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idAlbe, ncFRCid(idAlbe,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % albedoG)
      IF (exit_flag.ne.NoError) RETURN
# endif


# if (!defined ANA_CLOUD && !defined AIR_OCEAN) && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % cloudG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if (!defined ANA_SRFLUX && !defined AIR_OCEAN) && defined SHORTWAVE \
      && !defined ANA_NCEP
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % srflxG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
!
!-----------------------------------------------------------------------
!  Daily-averaged Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idAsrf, ncFRCid(idAsrf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % srflxG_avg)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined BULK_FLUXES && !defined LONGWAVE && !defined LONGWAVE_OUT \
   && !defined ANA_LRFLUX && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                      &
#  endif
     &                FORCES(ng) % lrflxG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined BULK_FLUXES && defined LONGWAVE_OUT &&                    \
    !defined ANA_LRFLUX && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % lrflxG)
      IF (exit_flag.ne.NoError) RETURN
# endif



# if (!defined ANA_RAIN && !defined AIR_OCEAN) && defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Rain fall rate.  Snow fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % rainG)
      IF (exit_flag.ne.NoError) RETURN
#  if defined SNOWFALL && !defined ANA_SNOW
      CALL get_2dfld (ng, iNLM, idsnow, ncFRCid(idsnow,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % snowG)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# if !defined ANA_RUNOFF && defined BULK_FLUXES && defined RUNOFF
!
!-----------------------------------------------------------------------
!  Surface runoff rate
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idRunoff, ncFRCid(idRunoff,ng),         &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % runoffG)
      IF (exit_flag.ne.NoError) RETURN
# endif


# if !defined ANA_STFLUX && !defined BULK_FLUXES && \
     !defined ATM2OCN_FLUXES && !defined NCEP_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % stflxG(:,:,:,itemp))
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
      defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface net heat flux from initial NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idTsur(itemp), BLK(ng)%ncid,          &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % stflxG(:,:,:,itemp))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# if !defined ANA_SST && defined QCORRECTION && !defined NCEP_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % sstG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_DQDSST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: heat flux sensitivity to
!  SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % dqdtG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % btflxG(:,:,:,itemp))
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SALINITY
#  if defined NL_BULK_FLUXES && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idEmPf, ncFRCid(idEmPf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % stflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN

#  elif !(defined ANA_SSFLUX || defined EMINUSP || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % stflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
       defined BULK_FLUXES && defined NL_BULK_FLUXES && \
       defined EMINUSP
!
!-----------------------------------------------------------------------
!  If not first NLM run, get freshwater flux (E-P) from initial NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idEmPf, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stflxG(:,:,:,isalt))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif

#  if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % sssG)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % btflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
#  ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stflxG(:,:,:,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  endif

#  ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTbot(i), ncFRCid(idTbot(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % btflxG(:,:,:,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  endif
# endif
#endif

#ifdef RED_TIDE
!
!-----------------------------------------------------------------------
!  Red tide Observed Dissolved Inorganic Nutrient.
!-----------------------------------------------------------------------
!
      CALL get_3dfld (ng, iNLM, idODIN, ncFRCid(idODIN,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                OCEAN(ng) % DIN_obsG)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.  In
!  grid refinement, only the coarser grid (RefineScale(ng)=0) open
!  boundary conditions data is processed and needed.
!=======================================================================
#ifndef ANA_FSOBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(iwest),                      &
     &                    ncBRYid(idZbry(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(ieast),                      &
     &                    ncBRYid(idZbry(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(isouth),                     &
     &                    ncBRYid(idZbry(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(inorth),                     &
     &                    ncBRYid(idZbry(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#endif

#ifndef ANA_M2OBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(iwest),                      &
     &                    ncBRYid(idU2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(iwest),                      &
     &                    ncBRYid(idV2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(ieast),                      &
     &                    ncBRYid(idU2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(ieast),                      &
     &                    ncBRYid(idV2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(isouth),                     &
     &                    ncBRYid(idU2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(isouth),                     &
     &                    ncBRYid(idV2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(inorth),                     &
     &                    ncBRYid(idU2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(inorth),                     &
     &                    ncBRYid(idV2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#endif

#ifdef SOLVE3D
# ifndef ANA_M3OBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest),                      &
     &                    ncBRYid(idU3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest),                      &
     &                    ncBRYid(idV3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast),                      &
     &                    ncBRYid(idU3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast),                      &
     &                    ncBRYid(idV3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth),                     &
     &                    ncBRYid(idU3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth),                     &
     &                    ncBRYid(idV3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth),                     &
     &                    ncBRYid(idU3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth),                     &
     &                    ncBRYid(idV3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifndef ANA_TOBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i),                  &
     &                      ncBRYid(idTbry(iwest,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i),                  &
     &                      ncBRYid(idTbry(ieast,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i),                 &
     &                      ncBRYid(idTbry(isouth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i),                 &
     &                      ncBRYid(idTbry(inorth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
# endif

# ifdef ICE_MODEL
#  ifndef ANA_MIOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idUibc(iwest),                      &
     &                    ncBRYid(idUibc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % uiG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(iwest,isVice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idVibc(iwest),                      &
     &                    ncBRYid(idVibc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % viG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isUice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idUibc(ieast),                      &
     &                    ncBRYid(idUibc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % uiG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isVice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idVibc(ieast),                      &
     &                    ncBRYid(idVibc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % viG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isUice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idUibc(isouth),                     &
     &                    ncBRYid(idUibc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % uiG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isVice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idVibc(isouth),                     &
     &                    ncBRYid(idVibc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % viG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isUice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idUibc(inorth),                     &
     &                    ncBRYid(idUibc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % uiG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isVice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idVibc(inorth),                     &
     &                    ncBRYid(idVibc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % viG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_AIOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isAice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idAibc(iwest),                      &
     &                    ncBRYid(idAibc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % aiG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isAice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idAibc(ieast),                      &
     &                    ncBRYid(idAibc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % aiG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isAice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idAibc(isouth),                     &
     &                    ncBRYid(idAibc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % aiG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isAice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idAibc(inorth),                     &
     &                    ncBRYid(idAibc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % aiG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif

#  ifndef ANA_HIOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isHice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHibc(iwest),                      &
     &                    ncBRYid(idHibc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hiG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isHice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHibc(ieast),                      &
     &                    ncBRYid(idHibc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hiG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isHice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHibc(isouth),                     &
     &                    ncBRYid(idHibc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hiG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isHice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHibc(inorth),                     &
     &                    ncBRYid(idHibc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hiG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_HSNOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isHsno,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHsnbc(iwest),                     &
     &                    ncBRYid(idHsnbc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hsnG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isHsno,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHsnbc(ieast),                     &
     &                    ncBRYid(idHsnbc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hsnG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isHsno,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHsnbc(isouth),                    &
     &                    ncBRYid(idHsnbc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hsnG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isHsno,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHsnbc(inorth),                    &
     &                    ncBRYid(idHsnbc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hsnG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_TIOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isTice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTibc(iwest),                      &
     &                    ncBRYid(idTibc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % tiG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isTice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTibc(ieast),                      &
     &                    ncBRYid(idTibc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % tiG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isTice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTibc(isouth),                     &
     &                    ncBRYid(idTibc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % tiG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isTice,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTibc(inorth),                     &
     &                    ncBRYid(idTibc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % tiG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  if defined MELT_PONDS && !defined ANA_SFWATOBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isApond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idApdbc(iwest),                     &
     &                    ncBRYid(idApdbc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % apondG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isApond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idApdbc(ieast),                     &
     &                    ncBRYid(idApdbc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % apondG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isApond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idApdbc(isouth),                    &
     &                    ncBRYid(idApdbc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % apondG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isApond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idApdbc(inorth),                    &
     &                    ncBRYid(idApdbc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % apondG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(iwest,isHpond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHpdbc(iwest),                     &
     &                    ncBRYid(idHpdbc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hpondG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isHpond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHpdbc(ieast),                     &
     &                    ncBRYid(idHpdbc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hpondG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isHpond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHpdbc(isouth),                    &
     &                    ncBRYid(idHpdbc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hpondG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isHpond,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idHpdbc(inorth),                    &
     &                    ncBRYid(idHpdbc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % hpondG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_SIG11OBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isSig11,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS11bc(iwest),                     &
     &                    ncBRYid(idS11bc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig11G_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isSig11,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS11bc(ieast),                     &
     &                    ncBRYid(idS11bc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig11G_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isSig11,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS11bc(isouth),                    &
     &                    ncBRYid(idS11bc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig11G_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isSig11,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS11bc(inorth),                    &
     &                    ncBRYid(idS11bc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig11G_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_SIG22OBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isSig22,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS22bc(iwest),                     &
     &                    ncBRYid(idS22bc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig22G_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isSig22,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS22bc(ieast),                     &
     &                    ncBRYid(idS22bc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig22G_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isSig22,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS22bc(isouth),                    &
     &                    ncBRYid(idS22bc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig22G_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isSig22,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS22bc(inorth),                    &
     &                    ncBRYid(idS22bc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig22G_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
#  ifndef ANA_SIG12OBC
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isSig12,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS12bc(iwest),                     &
     &                    ncBRYid(idS12bc(iwest),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig12G_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isSig12,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS12bc(ieast),                     &
     &                    ncBRYid(idS12bc(ieast),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig12G_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isSig12,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS12bc(isouth),                    &
     &                    ncBRYid(idS12bc(isouth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig12G_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isSig12,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idS12bc(inorth),                    &
     &                    ncBRYid(idS12bc(inorth),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % sig12G_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#  endif
# endif
#endif
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================

#ifndef ANA_SSH
!
!  Free-surface.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, ncCLMid(idSSHc,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  CLIMA(ng) % sshG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
#ifndef ANA_M2CLIMA
!
!  2D momentum.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, ncCLMid(idUbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  CLIMA(ng) % ubarclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, ncCLMid(idVbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  CLIMA(ng) % vbarclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
#ifdef SOLVE3D
# ifndef ANA_M3CLIMA
!
!  3D momentum.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, ncCLMid(idUclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  CLIMA(ng) % uclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, ncCLMid(idVclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  CLIMA(ng) % vclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif
# ifndef ANA_TCLIMA
!
!  Tracers.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iNLM, idTclm(i),                          &
     &                    ncCLMid(idTclm(i),ng),                        &
     &                    nCLMfiles(ng), CLM(1,ng), update(1),          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif

# if !defined ANA_AKTCLIMA && defined AKTCLIMATOLOGY
      CALL get_3dfld (ng, iNLM, idAclm, ncCLMid(idAclm,ng),             &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % AclmG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if (defined BEST_NPZ || defined BIO_GOANPZ) && defined IRON_LIMIT
      CALL get_3dfld (ng, iNLM, idTclm(iFe), ncCLMid(idTclm(iFe),ng),   &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % tclmG(:,:,:,:,iFe))
      IF (exit_flag.ne.NoError) RETURN
# endif
# if defined ICE_BIO && defined CLIM_ICE_1D
       CALL get_3dfld (ng, iNLM, idTclm(i1CI), ncCLMid(idTclm(i1CI),ng),&
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % tclmG(:,:,:,:,i1CI))
      IF (exit_flag.ne.NoError) RETURN
# endif
# if !defined ANA_OCLIMA && defined OCLIMATOLOGY
      CALL get_3dfld (ng, iNLM, idOclm, ncCLMid(idOclm,ng),             &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % oclmG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_AICLIMA && defined AICLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idAIclm, ncCLMid(idAIclm,ng),           &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % aiclmG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng, iNLM, idHIclm, ncCLMid(idHIclm,ng),           &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % hiclmG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_MICLIMA && defined MICLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idUIclm, ncCLMid(idUIclm,ng),           &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % uiclmG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng, iNLM, idVIclm, ncCLMid(idVIclm,ng),           &
     &                nCLMfiles(ng), CLM(1,ng), update(1),              &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                CLIMA(ng) % viclmG)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  If tangent linear model check, read in nonlinear forward solution
!  to compute dot product with perturbated nonlinear solution. Time
!  interpolation between snapshot is not required (see subroutine
!  "nl_dotproduct").
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
!
!  Read in free-surface.
!
        CALL get_2dfld (ng, iNLM, idFsur, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % zetaG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  OCEAN(ng) % ubarG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  OCEAN(ng) % vbarG)
        IF (exit_flag.ne.NoError) RETURN
# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % uG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), FWD(ng)%ncid,            &
     &                    1, FWD(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % tG(:,:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
      END IF
#endif

#if defined NLM_OUTER || defined TL_W4DPSAS           || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Read weak contraint forcing snapshots. Notice that the forward
!  basic state snapshops arrays are reused here.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse(ng)) THEN
!
!  Read in free-surface.
!
        CALL get_2dfld (ng, iNLM, idFsur, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % zetaG)
        IF (exit_flag.ne.NoError) RETURN
# ifndef SOLVE3D
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % ubarG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vbarG)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % uG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), TLF(ng)%ncid,            &
     &                    1, TLF(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % tG(:,:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
      END IF
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
#endif

      RETURN
      END SUBROUTINE get_data
