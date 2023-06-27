#include "cppdefs.h"
# ifdef TANGENT
      SUBROUTINE tl_get_data (ng)
!
!git $Id$
!svn $Id: tl_get_data.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
      logical, save :: recordless = .FALSE.

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
      CALL wclock_on (ng, iTLM, 3, __LINE__, MyFile)
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
        CALL get_ngfld (ng, iTLM, idRtra, SSF(ng)%ncid,                 &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                  SSF(ng)%pioFile,                                &
#  endif
     &                  1, SSF(ng), recordless, update(1),              &
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
          CALL get_ngfld (ng, iTLM, idRtrc(i), SSF(ng)%ncid,            &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    SSF(ng)%pioFile,                              &
#   endif
     &                    1, SSF(ng), recordless, update(1),            &
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
      CALL get_2dfld (ng, iTLM, idUair, FRCncid(idUair,ng),             &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idUair,ng),                            &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % UwindG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng , iTLM, idVair, FRCncid(idVair,ng),            &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idVair,ng),                            &
#  endif
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
      CALL get_2dfld (ng, iTLM, idUsms, FRCncid(idUsms,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idUsms,ng),                            &
#   endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                FORCES(ng) % sustrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idVsms, FRCncid(idVsms,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idVsms,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idUsms, BLK(ng)%ncid,                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#  endif
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FORCES(ng) % sustrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idVsms, BLK(ng)%ncid,                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#  endif
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                FORCES(ng) % svstrG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !defined ANA_PAIR       && \
     (defined ATM_PRESS      || defined BULK_FLUXES  || \
      defined ECOSIM)
#  if !(defined FRC_COUPLING || defined FORWARD_FLUXES)
!
!-----------------------------------------------------------------------
!  Surface air pressure from input FRC file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idPair, FRCncid(idPair,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idPair,ng),                            &
#   endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % PairG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  else
!
!-----------------------------------------------------------------------
!  Surface air pressure from NLM forward file.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idPair, BLK(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWdir, FRCncid(idWdir,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWdir,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWamp, FRCncid(idWamp,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWamp,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWlen, FRCncid(idWlen,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWlen,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWptp, FRCncid(idWptp,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWptp,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWpbt, FRCncid(idWpbt,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWpbt,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWorb, FRCncid(idWorb,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWorb,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWdis, FRCncid(idWdis,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWdis,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idWbrk, FRCncid(idWbrk,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idWbrk,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idCfra, FRCncid(idCfra,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idCfra,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idSrad, BLK(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#    endif
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
      CALL get_2dfld (ng, iTLM, idSrad, FRCncid(idSrad,ng),             &
#     if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idSrad,ng),                            &
#     endif
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
      CALL get_2dfld (ng, iTLM, idLrad, FRCncid(idLrad,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idLrad,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idLdwn, FRCncid(idLdwn,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idLdwn,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idTair, FRCncid(idTair,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idTair,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idQair, FRCncid(idQair,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idQair,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idrain, FRCncid(idrain,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idrain,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idTsur(itemp),                          &
     &                FRCncid(idTsur(itemp),ng),                        &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idTsur(itemp),ng),                     &
#    endif
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
      CALL get_2dfld (ng, iTLM, idTsur(itemp), BLK(ng)%ncid,            &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#    endif
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
      CALL get_2dfld (ng, iTLM, idSSTc, FRCncid(idSSTc,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idSSTc,ng),                            &
#   endif
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

      CALL get_2dfld (ng, iTLM, iddQdT, FRCncid(iddQdT,ng),             &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(iddQdT,ng),                            &
#   endif
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
      CALL get_2dfld (ng, iTLM, idTbot(itemp),                          &
     &                FRCncid(idTbot(itemp),ng),                        &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idTbot(itemp),ng),                     &
#   endif
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
      CALL get_2dfld (ng, iTLM, idsfwf, FRCncid(idsfwf,ng),             &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idsfwf,ng),                            &
#    endif
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
      CALL get_2dfld (ng, iTLM, idEmPf, BLK(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                BLK(ng)%pioFile,                                  &
#    endif
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
      CALL get_2dfld (ng, iTLM, idSSSc, FRCncid(idSSSc,ng),             &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idSSSc,ng),                            &
#    endif
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
      CALL get_2dfld (ng, iTLM, idTbot(isalt),                          &
     &                FRCncid(idTbot(isalt),ng),                        &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idTbot(isalt),ng),                     &
#    endif
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
        CALL get_2dfld (ng, iTLM, idTsur(i), FRCncid(idTsur(i),ng),     &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  FRCpioFile(idTsur(i),ng),                       &
#    endif
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
        CALL get_2dfld (ng, iTLM, idTbot(i), FRCncid(idTbot(i),ng),     &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  FRCpioFile(idTbot(i),ng),                       &
#    endif
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
          CALL get_ngfld (ng, iTLM, idZbry(iwest),                      &
     &                    BRYncid(idZbry(iwest),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idZbry(iwest),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idZbry(ieast),                      &
     &                    BRYncid(idZbry(ieast),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idZbry(ieast),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idZbry(isouth),                     &
     &                    BRYncid(idZbry(isouth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idZbry(isouth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idZbry(inorth),                     &
     &                    BRYncid(idZbry(inorth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idZbry(inorth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
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
          CALL get_ngfld (ng, iTLM, idU2bc(iwest),                      &
     &                    BRYncid(idU2bc(iwest),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU2bc(iwest),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV2bc(iwest),                      &
     &                    BRYncid(idV2bc(iwest),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV2bc(iwest),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU2bc(ieast),                      &
     &                    BRYncid(idU2bc(ieast),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU2bc(ieast),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV2bc(ieast),                      &
     &                    BRYncid(idV2bc(ieast),ng),                    &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV2bc(ieast),ng),                 &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU2bc(isouth),                     &
     &                    BRYncid(idU2bc(isouth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU2bc(isouth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV2bc(isouth),                     &
     &                    BRYncid(idV2bc(isouth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV2bc(isouth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU2bc(inorth),                     &
     &                    BRYncid(idU2bc(inorth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU2bc(inorth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV2bc(inorth),                     &
     &                    BRYncid(idV2bc(inorth),ng),                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV2bc(inorth),ng),                &
#  endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
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
          CALL get_ngfld (ng, iTLM, idU3bc(iwest),                      &
     &                    BRYncid(idU3bc(iwest),ng),                    &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU3bc(iwest),ng),                 &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV3bc(iwest),                      &
     &                    BRYncid(idV3bc(iwest),ng),                    &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV3bc(iwest),ng),                 &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU3bc(ieast),                      &
     &                    BRYncid(idU3bc(ieast),ng),                    &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU3bc(ieast),ng),                 &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV3bc(ieast),                      &
     &                    BRYncid(idV3bc(ieast),ng),                    &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV3bc(ieast),ng),                 &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU3bc(isouth),                     &
     &                    BRYncid(idU3bc(isouth),ng),                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU3bc(isouth),ng),                &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV3bc(isouth),                     &
     &                    BRYncid(idV3bc(isouth),ng),                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV3bc(isouth),ng),                &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idU3bc(inorth),                     &
     &                    BRYncid(idU3bc(inorth),ng),                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idU3bc(inorth),ng),                &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (tl_LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idV3bc(inorth),                     &
     &                    BRYncid(idV3bc(inorth),ng),                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    BRYpioFile(idV3bc(inorth),ng),                &
#   endif
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
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
            CALL get_ngfld (ng, iTLM, idTbry(iwest,i),                  &
     &                      BRYncid(idTbry(iwest,i),ng),                &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                      BRYpioFile(idTbry(iwest,i),ng),             &
#   endif
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (tl_LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iTLM, idTbry(ieast,i),                  &
     &                      BRYncid(idTbry(ieast,i),ng),                &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                      BRYpioFile(idTbry(ieast,i),ng),             &
#   endif
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (tl_LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iTLM, idTbry(isouth,i),                 &
     &                      BRYncid(idTbry(isouth,i),ng),               &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                      BRYpioFile(idTbry(isouth,i),ng),            &
#   endif
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (tl_LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iTLM, idTbry(inorth,i),                 &
     &                      BRYncid(idTbry(inorth,i),ng),               &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                      BRYpioFile(idTbry(inorth,i),ng),            &
#   endif
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
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
        CALL get_2dfld (ng, iTLM, idSSHc, CLMncid(idSSHc,ng),           &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                  CLMpioFile(idSSHc,ng),                          &
#  endif
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
        CALL get_2dfld (ng, iTLM, idUbcl, CLMncid(idUbcl,ng),           &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                  CLMpioFile(idUbcl,ng),                          &
#  endif
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  CLIMA(ng) % ubarclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iTLM, idVbcl, CLMncid(idVbcl,ng),           &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                  CLMpioFile(idVbcl,ng),                          &
#  endif
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
        CALL get_3dfld (ng, iTLM, idUclm, CLMncid(idUclm,ng),           &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                  CLMpioFile(idUclm,ng),                          &
#   endif
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % umask,                               &
#   endif
     &                  CLIMA(ng) % uclmG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_3dfld (ng, iTLM, idVclm, CLMncid(idVclm,ng),           &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                  CLMpioFile(idVclm,ng),                          &
#   endif
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
          CALL get_3dfld (ng, iTLM, idTclm(i),                          &
     &                    CLMncid(idTclm(i),ng),                        &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                    CLMpioFile(idTclm(i),ng),                     &
#   endif
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
      CALL get_2dfld (ng, iTLM, idFsur, FWD(ng)%ncid,                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#  endif
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
      CALL get_2dfld (ng, iTLM, idUbar, FWD(ng)%ncid,                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#  endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                OCEAN(ng) % ubarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL get_2dfld (ng, iTLM, idVbar, FWD(ng)%ncid,                   &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#  endif
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
      CALL get_2dfld (ng, iTLM, idRzet, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                OCEAN(ng) % rzetaG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idRu2d, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                OCEAN(ng) % rubarG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idRv2d, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
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
      CALL get_2dfld (ng, iTLM, idUfx1, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                COUPLING(ng) % DU_avg1G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idUfx2, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                COUPLING(ng) % DU_avg2G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idVfx1, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                COUPLING(ng) % DV_avg1G)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idVfx2, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
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
      CALL get_3dfld (ng, iTLM, idUvel, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                OCEAN(ng) % uG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_3dfld (ng, iRPM, idVvel, FWD(ng)%ncid,                   &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#   endif
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
      CALL get_2dfld (ng, iTLM, idRuct, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % umask,                                 &
#    endif
     &                COUPLING(ng) % rufrcG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_2dfld (ng, iTLM, idRvct, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#    endif
     &                COUPLING(ng) % rvfrcG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_3dfld (ng, iTLM, idRu3d, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % umask,                                 &
#    endif
     &                OCEAN(ng) % ruG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL get_3dfld (ng, iTLM, idRv3d, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
        CALL get_3dfld (ng, iTLM, idTvar(i), FWD(ng)%ncid,              &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                  FWD(ng)%pioFile,                                &
#   endif
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
        CALL get_3dfld (ng, iTLM, idDiff(i), FWD(ng)%ncid,              &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  FWD(ng)%pioFile,                                &
#    endif
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
     &                  MIXING(ng) % AktG(:,:,:,:,i))
        Fscale(idDiff(i),ng)=scale
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      scale=Fscale(idVvis,ng)                         ! save and rescale
      Fscale(idVvis,ng)=tl_Akv_fac(ng)
      CALL get_3dfld (ng, iTLM, idVvis, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_3dfld (ng, iTLM, idMtke, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_3dfld (ng, iTLM, idMtls, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_3dfld (ng, iTLM, idVmLS, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_3dfld (ng, iTLM, idVmKK, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_3dfld (ng, iTLM, idVmKP, FWD(ng)%ncid,                   &
#     if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#     endif
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
      CALL get_2dfld (ng, iTLM, idHsbl, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
      CALL get_2dfld (ng, iTLM, idHbbl, FWD(ng)%ncid,                   &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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
        CALL get_3dfld (ng, iTLM, idGhat(i), FWD(ng)%ncid,              &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                FWD(ng)%pioFile,                                  &
#    endif
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

#  ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Read frequent impulse forcing for weak constraint.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse(ng)) THEN
        CALL get_2dfld (ng, iTLM, idZtlf, TLF(ng)%ncid,                 &
#   if defined PIO_LIB && defined DISTRIBUTE
     &                  TLF(ng)%pioFile,                                &
#   endif
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
        CALL get_2dfld (ng, iTLM, idUbtf, TLF(ng)%ncid,                 &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  TLF(ng)%pioFile,                                &
#    endif
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % umask,                               &
#    endif
     &                  OCEAN(ng) % f_ubarG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_2dfld (ng, iTLM, idVbtf, TLF(ng)%ncid,                 &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  TLF(ng)%pioFile,                                &
#    endif
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
        CALL get_3dfld (ng, iTLM, idUtlf, TLF(ng)%ncid,                 &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  TLF(ng)%pioFile,                                &
#    endif
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % umask,                               &
#    endif
     &                  OCEAN(ng) % f_uG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL get_3dfld (ng, iTLM, idVtlf, TLF(ng)%ncid,                 &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                  TLF(ng)%pioFile,                                &
#    endif
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
          CALL get_3dfld (ng, iTLM, idTtlf(i), TLF(ng)%ncid,            &
#    if defined PIO_LIB && defined DISTRIBUTE
     &                    TLF(ng)%pioFile,                              &
#    endif
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
      CALL wclock_off (ng, iTLM, 3, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE tl_get_data
#else
      SUBROUTINE tl_get_data
      RETURN
      END SUBROUTINE tl_get_data
#endif
