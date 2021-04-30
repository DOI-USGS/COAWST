#include "cppdefs.h"
      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
#if defined HYPOXIA_SRM || defined RED_TIDE
      USE mod_biology
#endif
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined HYPOXIA_SRM              || \
    defined NLM_OUTER                || \
    defined RBL4DVAR                 || \
    defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined RED_TIDE                 || \
    defined SP4DVAR                  || \
    defined TLM_CHECK                || \
    defined TL_RBL4DVAR
      USE mod_ocean
#endif
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Lprocess

      logical, dimension(3) :: update = (/ .FALSE., .FALSE., .FALSE. /)
!
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
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
      CALL wclock_on (ng, iNLM, 3, __LINE__, MyFile)
#endif

#ifndef ANA_PSOURCE
!
!=======================================================================
!  Point Sources/Sinks time dependent data.
!=======================================================================
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
# endif
#endif
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file(s).
!=======================================================================
!
!  Set switch to process surface atmospheric fields.
!
#if defined FOUR_DVAR   && \
    defined BULK_FLUXES && defined PRIOR_BULK_FLUXES
!  In 4D-Var data assimilation applications, the user have the option
!  to fix the prior (background phase) surface fluxes in the successive
!  outer loops (Nouter>1) or the final analysis phase. In such case, the
!  fluxes are read from the background trajector
!
      IF (Nrun.eq.1) THEN
        Lprocess=.TRUE.
      ELSE
        Lprocess=.FALSE.
      END IF
#else
      Lprocess=.TRUE.
#endif

#if (!defined ANA_WINDS && !defined AIR_OCEAN) && \
    (defined BULK_FLUXES2D || defined BULK_FLUXES || \
     defined ECOSIM || defined SPECTRAL_LIGHT)
!
!  Surface wind components.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  FORCES(ng) % UwindG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),          &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  FORCES(ng) % VwindG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#if !(defined ANA_SMFLUX || defined FRC_COUPLING || defined BULK_FLUXES)
!
!  Surface wind stress components.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idUsms, ncFRCid(idUsms,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  FORCES(ng) % sustrG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iNLM, idVsms, ncFRCid(idVsms,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  FORCES(ng) % svstrG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#if !defined ANA_SSSFLX && defined SSSFLX
      CALL get_2dfld (ng, iNLM, idSSSf, ncFRCid(idSSSf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % sssflxG)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
#endif

#if !defined ANA_FASTICE && defined FASTICE_CLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idFastIce, ncFRCid(idFastIce,ng),       &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FORCES(ng) % fastice_clmG)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
#endif

#if (!defined ANA_PAIR && !defined AIR_OCEAN) && (defined BULK_FLUXES \
      || defined ECOSIM || defined ATM_PRESS || defined SPECTRAL_LIGHT)
!
!  Surface air pressure.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  FORCES(ng) % PairG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
#if !defined ANA_WWAVE && defined WAVE_DATA && !defined INWAVE_MODEL
!
!  Surface wind induced wave amplitude, direction and period.
!
# ifdef WAVES_DIR
      CALL get_2dfld (ng, iNLM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % DwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WAVES_DIRP
      CALL get_2dfld (ng, iNLM, idWdip, ncFRCid(idWdip,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % DwavepG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WAVES_HEIGHT
!
      CALL get_2dfld (ng, iNLM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % HwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WAVES_LENGTH
!
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % LwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef WAVES_LENGTHP
      CALL get_2dfld (ng, iNLM, idWlep, ncFRCid(idWlep,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % LwavepG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef WAVES_TOP_PERIOD
!
      CALL get_2dfld (ng, iNLM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_topG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WAVES_BOT_PERIOD
!
      CALL get_2dfld (ng, iNLM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_botG(:,:,1))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WAVES_UB
!
      CALL get_2dfld (ng, iNLM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Uwave_rmsG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined WAVES_DISS
      CALL get_2dfld (ng, iNLM, idWdib, ncFRCid(idWdib,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Dissip_breakG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iNLM, idWdiw, ncFRCid(idWdiw,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Dissip_wcapG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined ROLLER_SVENDSEN
      CALL get_2dfld (ng, iNLM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Wave_breakG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#ifdef SOLVE3D

# if !(defined ANA_CLOUD || defined FRC_COUPLING) && defined CLOUDS
!
!  Cloud fraction.
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % cloudG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !(defined ANA_SRFLUX || defined FRC_COUPLING) && defined SHORTWAVE
!
!  Surface solar shortwave radiation.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % srflxG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
!
!  Daily-averaged Surface solar shortwave radiation.
!
      CALL get_2dfld (ng, iNLM, idAsrf, ncFRCid(idAsrf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % srflxG_avg)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined BULK_FLUXES  && \
   !(defined FRC_COUPLING || defined LONGWAVE || defined LONGWAVE_OUT)
!
!  Surface net longwave radiation.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idLrad, ncFRCid(idLrad,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % lrflxG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if defined BULK_FLUXES  && defined LONGWAVE_OUT && \
    !defined FRC_COUPLING
!
!  Surface downwelling longwave radiation.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idLdwn, ncFRCid(idLdwn,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % lrflxG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if (!defined ANA_TAIR && !defined AIR_OCEAN) && \
    ( defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO) \
     || defined SPECTRAL_LIGHT)
!
!  Surface air temperature.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % TairG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if (!defined ANA_HUMIDITY && !defined AIR_OCEAN) && \
    ( defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO) \
    || defined SPECTRAL_LIGHT)
!
!  Surface air humidity.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % HairG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if !(defined ANA_RAIN || defined FRC_COUPLING) && defined BULK_FLUXES
!
!  Rain fall rate.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % rainG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!-----------------------------------------------------------------------
!  Snow fall rate.
!-----------------------------------------------------------------------
#  if defined SNOWFALL && !defined ANA_SNOW
      CALL get_2dfld (ng, iNLM, idsnow, ncFRCid(idsnow,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % snowG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

# endif

# if !(defined ANA_STFLUX || defined FRC_COUPLING || \
       defined BULK_FLUXES)
!
!  Surface net heat flux.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idTsur(itemp),                        &
     &                  ncFRCid(idTsur(itemp),ng),                      &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % stfluxG(:,:,:,itemp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if !defined ANA_SST && defined QCORRECTION
!
!  Surface net heat flux correction fields: sea surface temperature
!  (SST).
!
      CALL get_2dfld (ng, iNLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % sstG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !defined ANA_DQDSST && defined QCORRECTION
!
!  Surface net heat flux correction fields: heat flux sensitivity to
!  SST (dQdSST).
!
      CALL get_2dfld (ng, iNLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % dqdtG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifndef ANA_BTFLUX
!
!  Bottom net heat flux.
!
      CALL get_2dfld (ng, iNLM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % btfluxG(:,:,:,itemp))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !defined ANA_SSFLUX     && defined SALINITY
#  if !(defined BULK_FLUXES  || defined EMINUSP     || \
        defined FRC_COUPLING)
!
!  Surface net freshwater flux: E-P from NetCDF variable "swflux".
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idsfwf, ncFRCid(idsfwf,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stfluxG(:,:,:,isalt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  elif defined BULK_FLUXES  && !defined EMINUSP
!
!  Surface net freshwater flux (E-P) from NetCDF variable "EminusP".
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idEmPf, ncFRCid(idEMPf,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stfluxG(:,:,:,isalt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!  Surface net freshwater flux correction field: sea surface salinity.
!
      CALL get_2dfld (ng, iNLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % sssG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifndef ANA_BSFLUX
!
!  Bottom net freshwater flux.
!
      CALL get_2dfld (ng, iNLM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % btfluxG(:,:,:,isalt))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
#  ifndef ANA_SPFLUX
!
!  Passive tracers surface fluxes.
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stfluxG(:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#  endif

#  ifndef ANA_BPFLUX
!
!  Passive tracers bottom fluxes.
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTbot(i), ncFRCid(idTbot(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % btfluxG(:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#  endif
# endif
#endif

#if !defined ANA_RESPIRATION && defined HYPOXIA_SRM
!
!  Total respiration rate for hypoxia.
!
      CALL get_3dfld (ng, iNLM, idResR, ncFRCid(idResR,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                OCEAN(ng) % respirationG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef RED_TIDE
!
!  Red tide Observed Dissolved Inorganic Nutrient.
!
      CALL get_3dfld (ng, iNLM, idODIN, ncFRCid(idODIN,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                OCEAN(ng) % DIN_obsG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

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
!  The strategy is to save the 2D surface fluxes in the quicksave
!  NetCDF  file to allow frequent data records to resolve the daily
!  cycle while maintaining its size small.
!=======================================================================
!
      IF (.not.Lprocess) THEN
!
!  Get prior surface wind stress components.
!
        CALL get_2dfld (ng, iNLM, idUsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  FORCES(ng) % sustrG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iNLM, idVsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  FORCES(ng) % svstrG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ATM_PRESS
!
!  Get prior surface air pressure.
!
        CALL get_2dfld (ng, iNLM, idPair, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % PairG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  ifdef SHORTWAVE
!
!  Get prior shortwave radiation flux.  For consistency, we need to
!  process the same values using in the computation of the net heat
!  flux since there is a time interpolation from snapshots.
!
        CALL get_2dfld (ng, iNLM, idSrad, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % srflxG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Get prior surface net heat flux.
!
        CALL get_2dfld (ng, iNLM, idTsur(itemp), BLK(ng)%ncid,          &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % stfluxG(:,:,:,itemp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SALINITY
!
!  Get prior net freshwater flux (E-P). The kinematic net freshwater
!  will be computed in "set_vbc" by multiplying with surface salinity.
!
        CALL get_2dfld (ng, iNLM, idEmPf, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stfluxG(:,:,:,isalt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
      END IF
#endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================

#ifndef ANA_FSOBC
!
!  Free-surface open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(iwest),                      &
     &                    ncBRYid(idZbry(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(ieast),                      &
     &                    ncBRYid(idZbry(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(isouth),                     &
     &                    ncBRYid(idZbry(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(inorth),                     &
     &                    ncBRYid(idZbry(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
#endif

#ifndef ANA_M2OBC
!
!  2D momentum components open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(iwest),                      &
     &                    ncBRYid(idU2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(iwest),                      &
     &                    ncBRYid(idV2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(ieast),                      &
     &                    ncBRYid(idU2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(ieast),                      &
     &                    ncBRYid(idV2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(isouth),                     &
     &                    ncBRYid(idU2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(isouth),                     &
     &                    ncBRYid(idV2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(inorth),                     &
     &                    ncBRYid(idU2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(inorth),                     &
     &                    ncBRYid(idV2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
#endif

#ifdef SOLVE3D
# ifndef ANA_M3OBC
!
!  3D momentum components open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest),                      &
     &                    ncBRYid(idU3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest),                      &
     &                    ncBRYid(idV3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast),                      &
     &                    ncBRYid(idU3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast),                      &
     &                    ncBRYid(idV3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth),                     &
     &                    ncBRYid(idU3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth),                     &
     &                    ncBRYid(idV3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth),                     &
     &                    ncBRYid(idU3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth),                     &
     &                    ncBRYid(idV3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
# endif

# ifndef ANA_TOBC
!
!  Tracer variables open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i),                  &
     &                      ncBRYid(idTbry(iwest,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i),                  &
     &                      ncBRYid(idTbry(ieast,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i),                 &
     &                      ncBRYid(idTbry(isouth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i),                 &
     &                      ncBRYid(idTbry(inorth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
      END IF
# endif
#endif
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================

#ifndef ANA_SSH
!
!  Free-surface climatology.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, ncCLMid(idSSHc,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  CLIMA(ng) % sshG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
#ifndef ANA_M2CLIMA
!
!  2D momentum components climatology.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, ncCLMid(idUbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  CLIMA(ng) % ubarclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, ncCLMid(idVbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  CLIMA(ng) % vbarclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
#ifdef SOLVE3D
# ifndef ANA_M3CLIMA
!
!  3D momentum components climatology.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, ncCLMid(idUclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  CLIMA(ng) % uclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, ncCLMid(idVclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  CLIMA(ng) % vclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifndef ANA_TCLIMA
!
!  Tracers variables climatology.
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
# endif
#endif

#ifdef TLM_CHECK
!
!=======================================================================
!  If tangent linear model check, read in nonlinear forward solution
!  to compute dot product with perturbated nonlinear solution. Time
!  interpolation between snapshot is not required (see subroutine
!  "nl_dotproduct").
!=======================================================================
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  OCEAN(ng) % vbarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif
      END IF
#endif

#if defined NLM_OUTER                || \
    defined RBL4DVAR                 || \
    defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined SP4DVAR                  || \
    defined TL_RBL4DVAR
!
!=======================================================================
!  Read weak contraint forcing snapshots. Notice that the forward
!  basic state snapshops arrays are reused here.
!=======================================================================
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vbarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# else
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL wclock_off (ng, iNLM, 3, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE get_data
