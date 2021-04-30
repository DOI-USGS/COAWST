#include "cppdefs.h"
#ifdef TL_IOMS
      SUBROUTINE rp_get_data (ng)
!
!svn $Id: rp_get_data.F 1054 2021-03-06 19:47:12Z arango $
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
      USE mod_boundary
# ifdef FORWARD_READ
      USE mod_coupling
# endif
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
# ifdef FORWARD_READ
      USE mod_ocean
# endif
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
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
!
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile

# ifdef FORWARD_MIXING
!
      real(r8) :: scale
# endif
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

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iRPM, 3, __LINE__, MyFile)
# endif
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================

# ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iRPM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iRPM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), update(1),                        &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
# endif

# if !defined ANA_WINDS    && \
     ((defined BULK_FLUXES && !defined FORWARD_FLUXES) || \
      defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % UwindG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng , iRPM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % VwindG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifndef FRC_COUPLING
#  if !defined ANA_SMFLUX  && \
      !defined BULK_FLUXES && !defined FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components from input FRC file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                FORCES(ng) % sustrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                FORCES(ng) % svstrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined FORWARD_FLUXES || defined FRC_COUPLING
!
!-----------------------------------------------------------------------
!  Surface wind stress components from NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idUsms, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FORCES(ng) % sustrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idVsms, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                FORCES(ng) % svstrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !defined ANA_PAIR      && \
     (defined ECOSIM        || defined ATM_PRESS)
#  if !(defined BULK_FLUXES || defined FRC_COUPLING)
!
!-----------------------------------------------------------------------
!  Surface air pressure from input FRC file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % PairG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  elif defined BULK_FLUXES && defined FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Surface air pressure used in NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idPair, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % PairG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if !defined ANA_WWAVE && defined WAVE_DATA
!
!-----------------------------------------------------------------------
!  Surface wind induced wave properties.
!-----------------------------------------------------------------------
!
#  ifdef WAVES_DIR
      CALL get_2dfld (ng, iRPM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % DwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef WAVES_HEIGHT
!
      CALL get_2dfld (ng, iRPM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % HwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef WAVES_LENGTH
!
      CALL get_2dfld (ng, iRPM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % LwaveG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef WAVES_TOP_PERIOD
!
      CALL get_2dfld (ng, iRPM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % Pwave_topG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef WAVES_BOT_PERIOD
!
      CALL get_2dfld (ng, iRPM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % Pwave_botG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WAVES_UB
!
      CALL get_2dfld (ng, iRPM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % Ub_swanG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TKE_WAVEDISS
!
      CALL get_2dfld (ng, iRPM, idWdis, ncFRCid(idWdis,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % Wave_dissipG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SVENDSEN_ROLLER
!
      CALL get_2dfld (ng, iRPM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % Wave_breakG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef SOLVE3D

#  if !defined ANA_CLOUD && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % cloudG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SHORTWAVE && !defined ANA_SRFLUX
#   ifdef FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation from NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idSrad, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % srflxG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   else

#    if !(defined BULK_FLUXES || defined FRC_COUPLING)
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#     ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#     endif
     &                FORCES(ng) % srflxG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif

#   endif
#  endif

#  if (defined BULK_FLUXES && !defined FORWARD_FLUXES) && \
      !defined LONGWAVE    && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % lrflxG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if (defined BULK_FLUXES && !defined FORWARD_FLUXES) && \
      defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % lrflxG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if !defined ANA_TAIR && \
     ((defined BULK_FLUXES && !defined FORWARD_FLUXES) || \
      defined ECOSIM       || \
      (defined SHORTWAVE   && defined ANA_SRFLUX && defined ALBEDO))
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % TairG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if !defined ANA_HUMIDITY && \
      ((defined BULK_FLUXES && !defined FORWARD_FLUXES) || \
       defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % HairG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if !defined ANA_RAIN    && \
      (defined BULK_FLUXES && !defined FORWARD_FLUXES)
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % rainG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifndef ANA_STFLUX
#   if !(defined BULK_FLUXES || defined FORWARD_FLUXES)
!
!-----------------------------------------------------------------------
!  Surface net heat flux from input FRC file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % stfluxG(:,:,:,itemp))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   elif defined FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux from NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idTsur(itemp), BLK(ng)%ncid,            &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % stfluxG(:,:,:,itemp))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  if !defined ANA_SST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % sstG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if !defined ANA_DQDSST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: heat flux sensitivity to
!  SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % dqdtG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % btfluxG(:,:,:,itemp))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SALINITY        && !defined ANA_SSFLUX
#   if !(defined EMINUSP      || defined FORWARD_FLUXES || \
         defined FRC_COUPLING || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P from FRC file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % stfluxG(:,:,:,isalt))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   elif (defined EMINUSP      || defined FORWARD_FLUXES || \
          defined FRC_COUPLING)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idEmPf, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % stfluxG(:,:,:,isalt))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % sssG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % btfluxG(:,:,:,isalt))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
#   ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iRPM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
     &                  FORCES(ng) % stfluxG(:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#   endif

#   ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iRPM, idTbot(i), ncFRCid(idTbot(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
     &                  FORCES(ng) % btfluxG(:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#   endif
#  endif
# endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================

# ifndef ANA_FSOBC
!
      IF (LprocessOBC(ng)) THEN
        IF (tl_LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idZbry(iwest),                      &
     &                    ncBRYid(idZbry(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idZbry(ieast),                      &
     &                    ncBRYid(idZbry(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idZbry(isouth),                     &
     &                    ncBRYid(idZbry(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idZbry(inorth),                     &
     &                    ncBRYid(idZbry(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
# endif

# ifndef ANA_M2OBC
!
      IF (LprocessOBC(ng)) THEN
        IF (tl_LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU2bc(iwest),                      &
     &                    ncBRYid(idU2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV2bc(iwest),                      &
     &                    ncBRYid(idV2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU2bc(ieast),                      &
     &                    ncBRYid(idU2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV2bc(ieast),                      &
     &                    ncBRYid(idV2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU2bc(isouth),                     &
     &                    ncBRYid(idU2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV2bc(isouth),                     &
     &                    ncBRYid(idV2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU2bc(inorth),                     &
     &                    ncBRYid(idU2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV2bc(inorth),                     &
     &                    ncBRYid(idV2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
# endif

# ifdef SOLVE3D
#  ifndef ANA_M3OBC
!
      IF (LprocessOBC(ng)) THEN
        IF (tl_LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU3bc(iwest),                      &
     &                    ncBRYid(idU3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV3bc(iwest),                      &
     &                    ncBRYid(idV3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU3bc(ieast),                      &
     &                    ncBRYid(idU3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV3bc(ieast),                      &
     &                    ncBRYid(idV3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU3bc(isouth),                     &
     &                    ncBRYid(idU3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV3bc(isouth),                     &
     &                    ncBRYid(idV3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idU3bc(inorth),                     &
     &                    ncBRYid(idU3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        IF (tl_LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iRPM, idV3bc(inorth),                     &
     &                    ncBRYid(idV3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
#  endif

#  ifndef ANA_TOBC
!
      IF (LprocessOBC(ng)) THEN
        DO i=1,NT(ng)
          IF (tl_LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iRPM, idTbry(iwest,i),                  &
     &                      ncBRYid(idTbry(iwest,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (tl_LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iRPM, idTbry(ieast,i),                  &
     &                      ncBRYid(idTbry(ieast,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (tl_LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iRPM, idTbry(isouth,i),                 &
     &                      ncBRYid(idTbry(isouth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (tl_LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iRPM, idTbry(inorth,i),                 &
     &                      ncBRYid(idTbry(inorth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
      END IF
#  endif
# endif
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================

# ifndef ANA_SSH
!
!  Free-surface.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iRPM, idSSHc, ncCLMid(idSSHc,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  CLIMA(ng) % sshG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifndef ANA_M2CLIMA
!
!  2D momentum.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iRPM, idUbcl, ncCLMid(idUbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  CLIMA(ng) % ubarclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iRPM, idVbcl, ncCLMid(idVbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  CLIMA(ng) % vbarclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef SOLVE3D
#  ifndef ANA_M3CLIMA
!
!  3D momentum.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iRPM, idUclm, ncCLMid(idUclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % umask,                               &
#   endif
     &                  CLIMA(ng) % uclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_3dfld (ng, iRPM, idVclm, ncCLMid(idVclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#   endif
     &                  CLIMA(ng) % vclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifndef ANA_TCLIMA
!
!  Tracers.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iRPM, idTclm(i),                          &
     &                    ncCLMid(idTclm(i),ng),                        &
     &                    nCLMfiles(ng), CLM(1,ng), update(1),          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#   ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#   endif
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
# endif

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Read in forward state solution.
!-----------------------------------------------------------------------
!
!  Read in free-surface.
!
      CALL get_2dfld (ng, iRPM, idFsur, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                OCEAN(ng) % zetaG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read 2D momentum.
!
      CALL get_2dfld (ng, iRPM, idUbar, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                OCEAN(ng) % ubarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idVbar, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                OCEAN(ng) % vbarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FORWARD_RHS
!
!  Read in variables associated with 2D right-hand-side terms.
!
      CALL get_2dfld (ng, iRPM, idRzet, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                OCEAN(ng) % rzetaG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idRu2d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                OCEAN(ng) % rubarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idRv2d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                OCEAN(ng) % rvbarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
!
!  Read in variables associated with time-averaged 2D momentum terms.
!
      CALL get_2dfld (ng, iRPM, idUfx1, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                COUPLING(ng) % DU_avg1G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idUfx2, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                COUPLING(ng) % DU_avg2G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idVfx1, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                COUPLING(ng) % DV_avg1G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idVfx2, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                COUPLING(ng) % DV_avg2G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in 3D momentum.
!
      CALL get_3dfld (ng, iRPM, idUvel, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                OCEAN(ng) % uG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_3dfld (ng, iRPM, idVvel, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                OCEAN(ng) % vG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef FORWARD_RHS
!
!  Read in variables associated with 3D momentum right-hand-side terms.
!
      CALL get_2dfld (ng, iRPM, idRuct, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % umask,                                 &
#    endif
     &                COUPLING(ng) % rufrcG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iRPM, idRvct, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#    endif
     &                COUPLING(ng) % rvfrcG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_3dfld (ng, iRPM, idRu3d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % umask,                                 &
#    endif
     &                OCEAN(ng) % ruG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_3dfld (ng, iRPM, idRv3d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#    endif
     &                OCEAN(ng) % rvG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
!
!  Read in 3D tracers.
!
      DO i=1,NT(ng)
        CALL get_3dfld (ng, iRPM, idTvar(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  OCEAN(ng) % tG(:,:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

#   ifdef FORWARD_MIXING
!
!  Read in vertical mixing variables.
!
      DO i=1,NAT
        scale=Fscale(idDiff(i),ng)                    ! save and rescale
        Fscale(idDiff(i),ng)=tl_Akt_fac(i,ng)
        CALL get_3dfld (ng, iRPM, idDiff(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
     &                  MIXING(ng) % AktG(:,:,:,:,i))
        Fscale(idDiff(i),ng)=scale
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

      scale=Fscale(idVvis,ng)                         ! save and rescale
      Fscale(idVvis,ng)=tl_Akv_fac(ng)
      CALL get_3dfld (ng, iRPM, idVvis, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % AkvG)
      Fscale(idVvis,ng)=scale
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   if defined MY25_MIXING_NOT_YET || defined GLS_MIXING_NOT_YET
!
!  Read in turbulent kinetic energy.
!
      CALL get_3dfld (ng, iRPM, idMtke, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % tkeG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in turbulent kinetic energy times length scale.
!
      CALL get_3dfld (ng, iRPM, idMtls, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % glsG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in vertical mixing length scale.
!
      CALL get_3dfld (ng, iRPM, idVmLS, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % LscaleG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in vertical mixing coefficient for turbulent kinetic energy.
!
      CALL get_3dfld (ng, iRPM, idVmKK, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % AkkG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef GLS_MIXING_NOT_YET
!
!  Read in vertical mixing coefficient for turbulent length scale.
!
      CALL get_3dfld (ng, iRPM, idVmKP, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#     ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#     endif
     &                MIXING(ng) % AkpG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif

#   ifdef LMD_MIXING_NOT_YET
!
!  Read in depth of surface oceanic boundary layer.
!
      CALL get_2dfld (ng, iRPM, idHsbl, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % hsblG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef LMD_BKPP_NOT_YET
!
!  Read in depth of bottom oceanic boundary layer.
!
      CALL get_2dfld (ng, iRPM, idHbbl, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                MIXING(ng) % hbblG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef LMD_NONLOCAL_NOT_YET
!
!  Read in boundary layer nonlocal transport.
!
      DO i=1,NAT
        CALL get_3dfld (ng, iRPM, idGhat(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
     &                  MIXING(ng) % ghatsG(:,:,:,i))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#   endif
#  endif

#  if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
      defined TL_R4DVAR
!
!-----------------------------------------------------------------------
!  Read frequent impulse forcing for weak constraint.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse(ng)) THEN
        CALL get_2dfld (ng, iRPM, idZtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  OCEAN(ng) % f_zetaG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifndef SOLVE3D
!
!  Read 2D momentum forcing.
!
        CALL get_2dfld (ng, iRPM, idUbtf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % umask,                               &
#    endif
     &                  OCEAN(ng) % f_ubarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_2dfld (ng, iRPM, idVbtf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#    endif
     &                  OCEAN(ng) % f_vbarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef SOLVE3D
!
!  Read in 3D momentum forcing.
!
        CALL get_3dfld (ng, iRPM, idUtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % umask,                               &
#    endif
     &                  OCEAN(ng) % f_uG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_3dfld (ng, iRPM, idVtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#    endif
     &                  OCEAN(ng) % f_vG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in 3D tracers forcing.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iRPM, idTtlf(i), TLF(ng)%ncid,            &
     &                    1, TLF(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#    ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#    endif
     &                    OCEAN(ng) % f_tG(:,:,:,:,i))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
#   endif
      END IF
#  endif
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iRPM, 3, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE rp_get_data
#else
      SUBROUTINE rp_get_data
      RETURN
      END SUBROUTINE rp_get_data
#endif
