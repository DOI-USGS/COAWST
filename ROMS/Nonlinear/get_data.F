#include "cppdefs.h"
      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 800 2008-10-21 00:17:29Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and assimilation data    !
!  from input NetCDF files.  If there is more than one time-record,    !
!  data  is loaded  into global two-time record arrays.  The actual    !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
#ifdef CLIMATOLOGY
      USE mod_clima
#endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined ASSIMILATION || defined NUDGING
      USE mod_obs
#endif
#if defined NLM_OUTER || defined TLM_CHECK || defined TL_W4DPSAS || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
      USE mod_ocean
#endif
      USE mod_scalars
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources
#endif
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

#ifdef OBC
      integer :: ILB, IUB, JLB, JUB
#endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i

#ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
#endif
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
#ifdef FRC_FILE
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
#endif
#if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                             defined TS_PSOURCE || defined Q_PSOURCE)
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
# if defined UV_PSOURCE || defined Q_PSOURCE
!
!  Point Source/Sink vertically integrated mass transport.
!
      CALL get_ngfld (ng, iNLM, idRtra, ncFRCid(idRtra,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,                &
     &                SOURCES(ng) % QbarG(1,1))
# endif
# if defined TS_PSOURCE && defined SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), ncFRCid(idRtrc(i),ng),   &
     &                    nFfiles(ng), FRCname(1,ng), update(1),        &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(1,1,1,i))
        END IF
      END DO
# endif
#endif
#if (!defined ANA_WINDS && !defined AIR_OCEAN) && \
    (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                FORCES(ng) % UwindG(LBi,LBj,1))
      CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                FORCES(ng) % VwindG(LBi,LBj,1))
#endif

#ifndef AIR_OCEAN
# if !defined ANA_SMFLUX && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % sustrG(LBi,LBj,1))
      CALL get_2dfld (ng, iNLM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % svstrG(LBi,LBj,1))
# endif

# if defined FOUR_DVAR && defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface wind stress components from initial
!  NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idUsms, ncBLKid(ng), 1,               &
     &                  BLKname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#  endif
     &                  FORCES(ng) % sustrG(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idVsms, ncBLKid(ng), 1,               &
     &                  BLKname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#  endif
     &                  FORCES(ng) % svstrG(LBi,LBj,1))
      END IF
# endif
#endif
#if !defined ANA_WWAVE && defined WAVE_DATA && !defined INWAVE_COUPLING
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
# ifdef WAVES_DIR
      CALL get_2dfld (ng, iNLM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % DwaveG(LBi,LBj,1))
# endif
# ifdef WAVES_HEIGHT
      CALL get_2dfld (ng, iNLM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % HwaveG(LBi,LBj,1))
# endif
# ifdef WAVES_LENGTH
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % LwaveG(LBi,LBj,1))
# endif
# ifdef WAVES_LENGTHP
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlep,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % LwavepG(LBi,LBj,1))
# endif
# ifdef WAVES_TOP_PERIOD
      CALL get_2dfld (ng, iNLM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Pwave_topG(LBi,LBj,1))
# endif
# ifdef WAVES_BOT_PERIOD
      CALL get_2dfld (ng, iNLM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Pwave_botG(LBi,LBj,1))
# endif
# if defined WAVES_UB
      CALL get_2dfld (ng, iNLM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Uwave_rmsG(LBi,LBj,1))
# endif
# if defined WAVES_DISS
      CALL get_2dfld (ng, iNLM, idWdib, ncFRCid(idWdib,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Dissip_breakG(LBi,LBj,1))
      CALL get_2dfld (ng, iNLM, idWdiw, ncFRCid(idWdiw,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Dissip_wcapG(LBi,LBj,1))
# endif
# if defined ROLLER_SVENDSEN
      CALL get_2dfld (ng, iNLM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % Wave_breakG(LBi,LBj,1))
# endif
#endif

#ifdef SOLVE3D

# if (!defined ANA_CLOUD && !defined AIR_OCEAN) && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % cloudG(LBi,LBj,1))
# endif

# if (!defined ANA_SRFLUX && !defined AIR_OCEAN) && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % srflxG(LBi,LBj,1))
# endif

# if defined BULK_FLUXES && !defined LONGWAVE && !defined LONGWAVE_OUT \
   && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
# endif

# if defined BULK_FLUXES && defined LONGWAVE_OUT && !defined AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
# endif

# if (!defined ANA_PAIR && !defined AIR_OCEAN) && (defined BULK_FLUXES \
      || defined ECOSIM || defined ATM_PRESS)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % PairG(LBi,LBj,1))
# endif

# if (!defined ANA_TAIR && !defined AIR_OCEAN) && \
     ( defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO) )
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % TairG(LBi,LBj,1))
# endif

# if (!defined ANA_HUMIDITY && !defined AIR_OCEAN) && \
     (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % HairG(LBi,LBj,1))
# endif

# if (!defined ANA_RAIN && !defined AIR_OCEAN) && defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % rainG(LBi,LBj,1))
# endif

# if !defined ANA_STFLUX && !defined BULK_FLUXES && \
     !defined ATM2OCN_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,itemp))
# endif

# if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
      defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface net heat flux from initial NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idTsur(itemp), ncBLKid(ng), 1,        &
     &                  BLKname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  FORCES(ng) % stflxG(LBi,LBj,1,itemp))
      END IF
# endif

# if !defined ANA_SST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST) and heat flux sensitivity to SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % sstG(LBi,LBj,1))
      CALL get_2dfld (ng, iNLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % dqdtG(LBi,LBj,1))
# endif

# ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,itemp))
# endif

# ifdef SALINITY
#  if defined NL_BULK_FLUXES && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idEmPf, ncFRCid(idEmPf,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,isalt))

#  elif !(defined ANA_SSFLUX || defined EMINUSP || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,isalt))
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
        CALL get_2dfld (ng, iNLM, idEmPf, ncBLKid(ng), 1,               &
     &                  BLKname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  FORCES(ng) % stflxG(LBi,LBj,1,isalt))
      END IF
#  endif

#  if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % sssG(LBi,LBj,1))
#  endif

#  ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,isalt))
#  endif
# endif

# if defined SEDIMENT || defined BIOLOGY
#  ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  FORCES(ng) % stflxG(LBi,LBj,1,i))
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
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  FORCES(ng) % btflxG(LBi,LBj,1,i))
      END DO
#  endif
# endif
#endif

#ifdef OBC_DATA
# ifdef REFINED_GRID
      IF (ng.eq.1) THEN
# endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================
!
# ifndef ANA_FSOBC
#  ifdef WEST_FSOBC
      CALL get_ngfld (ng, iNLM, idZbry(iwest), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % zetaG_west(JLB,1))
#  endif
#  ifdef EAST_FSOBC
      CALL get_ngfld (ng, iNLM, idZbry(ieast), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % zetaG_east(JLB,1))
#  endif
#  ifdef SOUTH_FSOBC
      CALL get_ngfld (ng, iNLM, idZbry(isouth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % zetaG_south(ILB,1))
#  endif
#  ifdef NORTH_FSOBC
      CALL get_ngfld (ng, iNLM, idZbry(inorth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % zetaG_north(ILB,1))
#  endif
# endif
# ifndef ANA_M2OBC
#  ifdef WEST_M2OBC
      CALL get_ngfld (ng, iNLM, idU2bc(iwest), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % ubarG_west(JLB,1))
      CALL get_ngfld (ng, iNLM, idV2bc(iwest), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % vbarG_west(JLB,1))
#  endif
#  ifdef EAST_M2OBC
      CALL get_ngfld (ng, iNLM, idU2bc(ieast), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % ubarG_east(JLB,1))
      CALL get_ngfld (ng, iNLM, idV2bc(ieast), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % vbarG_east(JLB,1))
#  endif
#  ifdef SOUTH_M2OBC
      CALL get_ngfld (ng, iNLM, idU2bc(isouth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % ubarG_south(ILB,1))
      CALL get_ngfld (ng, iNLM, idV2bc(isouth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % vbarG_south(ILB,1))
#  endif
#  ifdef NORTH_M2OBC
      CALL get_ngfld (ng, iNLM, idU2bc(inorth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % ubarG_north(ILB,1))
      CALL get_ngfld (ng, iNLM, idV2bc(inorth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                   &
     &                BOUNDARY(ng) % vbarG_north(ILB,1))
#  endif
# endif
# ifdef SOLVE3D
#  ifndef ANA_M3OBC
#   ifdef WEST_M3OBC
      CALL get_ngfld (ng, iNLM, idU3bc(iwest), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % uG_west(JLB,1,1))
      CALL get_ngfld (ng, iNLM, idV3bc(iwest), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % vG_west(JLB,1,1))
#   endif
#   ifdef EAST_M3OBC
      CALL get_ngfld (ng, iNLM, idU3bc(ieast), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % uG_east(JLB,1,1))
      CALL get_ngfld (ng, iNLM, idV3bc(ieast), ncBRYid(ng), 1,          &
     &                BRYname(ng), update(1),                           &
     &                JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % vG_east(JLB,1,1))
#   endif
#   ifdef SOUTH_M3OBC
      CALL get_ngfld (ng, iNLM, idU3bc(isouth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % uG_south(ILB,1,1))
      CALL get_ngfld (ng, iNLM, idV3bc(isouth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % vG_south(ILB,1,1))
#   endif
#   ifdef NORTH_M3OBC
      CALL get_ngfld (ng, iNLM, idU3bc(inorth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % uG_north(ILB,1,1))
      CALL get_ngfld (ng, iNLM, idV3bc(inorth), ncBRYid(ng), 1,         &
     &                BRYname(ng), update(1),                           &
     &                ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),           &
     &                BOUNDARY(ng) % vG_north(ILB,1,1))
#   endif
#  endif
#  ifndef ANA_TOBC
#   ifdef WEST_TOBC
#    ifdef NOSED_OBC
        DO i=1,NAT
#    else
        DO i=1,NT(ng)
#    endif
        CALL get_ngfld (ng, iNLM, idTbry(iwest,i), ncBRYid(ng), 1,      &
     &                  BRYname(ng), update(1),                         &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % tG_west(JLB,1,1,i))
      END DO
#   endif
#   ifdef EAST_TOBC
#    ifdef NOSED_OBC
        DO i=1,NAT
#    else
        DO i=1,NT(ng)
#    endif
        CALL get_ngfld (ng, iNLM, idTbry(ieast,i), ncBRYid(ng), 1,      &
     &                  BRYname(ng), update(1),                         &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % tG_east(JLB,1,1,i))
      END DO
#   endif
#   ifdef SOUTH_TOBC
#    ifdef NOSED_OBC
        DO i=1,NAT
#    else
        DO i=1,NT(ng)
#    endif
        CALL get_ngfld (ng, iNLM, idTbry(isouth,i), ncBRYid(ng), 1,     &
     &                  BRYname(ng), update(1),                         &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % tG_south(ILB,1,1,i))
      END DO
#   endif
#   ifdef NORTH_TOBC
#    ifdef NOSED_OBC
        DO i=1,NAT
#    else
        DO i=1,NT(ng)
#    endif
        CALL get_ngfld (ng, iNLM, idTbry(inorth,i), ncBRYid(ng), 1,     &
     &                  BRYname(ng), update(1),                         &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % tG_north(ILB,1,1,i))
      END DO
#   endif
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
#endif

#ifdef CLM_FILE
!
!=======================================================================
!  Read in climatology data from CLIMATOLOGY NetCDF file.
!=======================================================================
!
# if !defined ANA_SSH && defined ZCLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idSSHc, ncCLMid(ng), 1,                 &
     &                CLMname(ng), update(1),                           &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % sshG(LBi,LBj,1))
# endif
# if !defined ANA_M2CLIMA && defined M2CLIMATOLOGY
      CALL get_2dfld (ng, iNLM, idUbcl, ncCLMid(ng), 1,                 &
     &                CLMname(ng), update(1),                           &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % ubarclmG(LBi,LBj,1))
      CALL get_2dfld (ng, iNLM, idVbcl, ncCLMid(ng), 1,                 &
     &                CLMname(ng), update(1),                           &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % vbarclmG(LBi,LBj,1))
# endif
# ifdef SOLVE3D
#  if !defined ANA_TCLIMA && defined TCLIMATOLOGY
      DO i=1,NAT
        CALL get_3dfld (ng, iNLM, idTclm(i), ncCLMid(ng), 1,            &
     &                  CLMname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  CLIMA(ng) % tclmG(LBi,LBj,1,1,i))
      END DO
#  endif
#  if !defined ANA_M3CLIMA && defined M3CLIMATOLOGY
      CALL get_3dfld (ng, iNLM, idUclm, ncCLMid(ng), 1,                 &
     &                CLMname(ng), update(1),                           &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % uclmG(LBi,LBj,1,1))
      CALL get_3dfld (ng, iNLM, idVclm, ncCLMid(ng), 1,                 &
     &                CLMname(ng), update(1),                           &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % vclmG(LBi,LBj,1,1))
#  endif
# endif
#endif

#if defined ASSIMILATION || defined NUDGING
!
!=======================================================================
!  Read in assimilation data from ASSIMILATION NetCDF files.
!=======================================================================
!
# ifdef ASSIMILATION_SSH
!
!-----------------------------------------------------------------------
!  Read in sea surface height observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SSH(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHo, ncSSHid(ng), 1,               &
     &                  SSHname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  OBS(ng) % SSHobs(LBi,LBj))
        CALL get_2dfld (ng, iNLM, idSSHe, ncSSHid(ng), 1,               &
     &                  SSHname(ng), update(2),                         &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  OBS(ng) % EobsSSH(LBi,LBj))
        IF (update(1).and.update(2)) THEN
          update_SSH(ng)=.TRUE.
          tSSHobs(1,ng)=Vtime(1,idSSHo,ng)
          tsSSHobs(ng)=Vtime(1,idSSHo,ng)*day2sec
          Finfo(7,idSSHo,ng)=tsSSHobs(ng)
          Finfo(7,idSSHe,ng)=tsSSHobs(ng)
          EobsSSHmin(ng)=Finfo(8,idSSHe,ng)
          EobsSSHmax(ng)=Finfo(9,idSSHe,ng)
        END IF
      END IF
# endif

# if defined NUDGING_SST || defined ASSIMILATION_SST
!
!-----------------------------------------------------------------------
!  Read in sea surface temperature observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SST(ng)) THEN
#  ifdef NUDGING_SST
        CALL get_2dfld (ng, iNLM, idSSTo, ncSSTid(ng), 1,               &
     &                  SSTname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSTdat(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idSSTe, ncSSTid(ng), 1,               &
     &                  SSTname(ng), update(2),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatSST(LBi,LBj,1))
#  else
        CALL get_2dfld (ng, iNLM, idSSTo, ncSSTid(ng), 1,               &
     &                  SSTname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSTobs(LBi,LBj))
        CALL get_2dfld (ng, iNLM, idSSTe, ncSSTid(ng), 1,               &
     &                  SSTname(ng), update(2),                         &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsSST(LBi,LBj))
#  endif
        IF (update(1).and.update(2)) update_SST(ng)=.TRUE.
      END IF
# endif

# if defined NUDGING_T || defined ASSIMILATION_T
!
!-----------------------------------------------------------------------
!  Read in tracers observations and error variance.
!-----------------------------------------------------------------------
!
      DO i=1,NAT
        IF (assi_T(i,ng)) THEN
#  ifdef NUDGING_T
          CALL get_3dfld (ng, iNLM, idTobs(i), ncTOBSid(ng), 1,         &
     &                    TOBSname(ng), update(1),                      &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % Tdat(LBi,LBj,1,1,i))
          CALL get_3dfld (ng, iNLM, idTerr(i), ncTOBSid(ng), 1,         &
     &                    TOBSname(ng), update(2),                      &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % EdatT(LBi,LBj,1,1,i))
          IF (update(1).and.update(2)) update_T(i,ng)=.TRUE.
#  else
          CALL get_3dfld (ng, iNLM, idTobs(i), ncTOBSid(ng), 1,         &
     &                    TOBSname(ng), update(1),                      &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,           &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % Tobs(LBi,LBj,1,i))
          CALL get_3dfld (ng, iNLM, idTerr(i), ncTOBSid(ng), 1,         &
     &                    TOBSname(ng), update(2),                      &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,           &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % EobsT(LBi,LBj,1,i))
          IF (update(1).and.update(2)) THEN
            update_T(i,ng)=.TRUE.
            tTobs(1,i,ng)=Vtime(1,idTobs(i),ng)
            tsTobs(i,ng)=Vtime(1,idTobs(i),ng)*day2sec
            Finfo(7,idTobs(i),ng)=tsTobs(i,ng)
            Finfo(7,idTerr(i),ng)=tsTobs(i,ng)
            EobsTmin(i,ng)=Finfo(8,idTerr(i),ng)
            EobsTmax(i,ng)=Finfo(9,idTerr(i),ng)
          END IF
#  endif
        END IF
      END DO
# endif

# if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
!
!-----------------------------------------------------------------------
!  Read in surface current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UVsur(ng)) THEN
#  ifdef NUDGING_UVsur
        CALL get_2dfld (ng, iNLM, idUsur, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(1),                        &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Usurdat(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idVsur, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(2),                        &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vsurdat(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idUVse, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(3),                        &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatVsur(LBi,LBj,1))
#  else
        CALL get_2dfld (ng, iNLM, idUsur, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(1),                        &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Usur(LBi,LBj))
        CALL get_2dfld (ng, iNLM, idVsur, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(2),                        &
     &                  LBi, UBi, LBj, UBj, 1 ,1,                       &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vsur(LBi,LBj))
        CALL get_2dfld (ng, iNLM, idUVse, ncVSURid(ng), 1,              &
     &                  VSURname(ng), update(3),                        &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsVsur(LBi,LBj))
#  endif
        IF (update(1).and.update(2).and.update(3)) update_UVsur=.TRUE.
      END IF
# endif

# if defined NUDGING_UV || defined ASSIMILATION_UV
!
!-----------------------------------------------------------------------
!  Read in horizontal current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UV(ng)) THEN
#  ifdef NUDGING_UV
        CALL get_3dfld (ng, iNLM, idUobs, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(1),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Udat(LBi,LBj,1,1))
        CALL get_3dfld (ng, iNLM, idVobs, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(2),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vdat(LBi,LBj,1,1))
        CALL get_3dfld (ng, iNLM, idUVer, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(3),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatUV(LBi,LBj,1,1))
        IF (update(1).and.update(2).and.update(3)) update_UV(ng)=.TRUE.
#  else
        CALL get_3dfld (ng, iNLM, idUobs, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(1),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Uobs(LBi,LBj,1))
        CALL get_3dfld (ng, iNLM, idVobs, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(2),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vobs(LBi,LBj,1))
        CALL get_3dfld (ng, iNLM, idUVer, ncVOBSid(ng), 1,              &
     &                  VOBSname(ng), update(3),                        &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsUV(LBi,LBj,1))
        IF (update(1).and.update(2).and.update(3)) THEN
          update_UV(ng)=.TRUE.
          tVobs(1,ng)=Vtime(1,idVobs,ng)
          tsVobs(ng)=Vtime(1,idVobs,ng)*day2sec
          Finfo(7,idUobs,ng)=tsVobs(ng)
          Finfo(7,idVobs,ng)=tsVobs(ng)
          Finfo(7,idUVer,ng)=tsVobs(ng)
          EobsUVmin=Finfo(8,idUVer,ng)
          EobsUVmax=Finfo(9,idUVer,ng)
        END IF
#  endif
      END IF
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
        CALL get_2dfld (ng, iNLM, idFsur, ncFWDid(ng), 1,               &
     &                  FWDname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % zetaG(LBi,LBj,1))
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, ncFWDid(ng), 1,               &
     &                  FWDname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % ubarG(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idVbar, ncFWDid(ng), 1,               &
     &                  FWDname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % vbarG(LBi,LBj,1))
# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, ncFWDid(ng), 1,               &
     &                  FWDname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % uG(LBi,LBj,1,1))
        CALL get_3dfld (ng, iNLM, idVvel, ncFWDid(ng), 1,               &
     &                  FWDname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % vG(LBi,LBj,1,1))
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), ncFWDid(ng), 1,          &
     &                    FWDname(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % tG(LBi,LBj,1,1,i))

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
      IF (FrequentImpulse) THEN
!
!  Read in free-surface.
!
        CALL get_2dfld (ng, iNLM, idFsur, ncTLFid(ng), 1,               &
     &                  TLFname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % zetaG(LBi,LBj,1))
# ifndef SOLVE3D
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, ncTLFid(ng), 1,               &
     &                  TLFname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % ubarG(LBi,LBj,1))
        CALL get_2dfld (ng, iNLM, idVbar, ncTLFid(ng), 1,               &
     &                  TLFname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % vbarG(LBi,LBj,1))
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, ncTLFid(ng), 1,               &
     &                  TLFname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % uG(LBi,LBj,1,1))
        CALL get_3dfld (ng, iNLM, idVvel, ncTLFid(ng), 1,               &
     &                  TLFname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#  endif
     &                  OCEAN(ng) % vG(LBi,LBj,1,1))
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), ncTLFid(ng), 1,          &
     &                    TLFname(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % tG(LBi,LBj,1,1,i))
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
