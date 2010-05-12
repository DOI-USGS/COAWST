#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_get_data (ng)
!
!svn $Id: ad_get_data.F 451 2010-03-04 21:43:03Z arango $
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
# if defined FORWARD_READ && defined SOLVE3D
      USE mod_coupling
# endif
# if defined AD_SENSITIVITY      || defined CLIMATOLOGY      || \
     defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SENSITIVITY_4DVAR   || defined SO_SEMI
      USE mod_clima
# endif
      USE mod_forces
# ifdef SENSITIVITY_4DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
# if defined ASSIMILATION || defined NUDGING
      USE mod_obs
# endif
# ifdef FORWARD_READ
      USE mod_ocean
# endif
      USE mod_scalars
# if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
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

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i

# ifdef FORWARD_MIXING
      real(r8) :: scale
# endif

# ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif
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
      CALL wclock_on (ng, iADM, 3)
# endif
# ifdef FRC_FILE
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
# endif

# if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                              defined TS_PSOURCE || defined Q_PSOURCE)
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile.  Point sources are at U-
!  and V-points.
!-----------------------------------------------------------------------
!
#  if defined UV_PSOURCE || defined Q_PSOURCE
!
!  Point Source/Sink vertically integrated mass transport.
!
      CALL get_ngfldr (ng, iADM, idRtra, ncFRCid(idRtra,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,               &
     &                 SOURCES(ng) % QbarG(1,1))
#  endif
#  if defined TS_PSOURCE && defined SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfldr (ng, iADM, idRtrc(i), ncFRCid(idRtrc(i),ng),  &
     &                     nFfiles(ng), FRCname(1,ng), update(1),       &
     &                     1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),   &
     &                     SOURCES(ng) % TsrcG(1,1,1,i))
        END IF
      END DO
#  endif
# endif

# if !defined ANA_WINDS    && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idUair, ncFRCid(idUair,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 FORCES(ng) % UwindG(LBi,LBj,1))
      CALL get_2dfldr (ng , iADM, idVair, ncFRCid(idVair,ng),           &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 FORCES(ng) % VwindG(LBi,LBj,1))
# endif

# ifndef AIR_OCEAN
#  if !defined ANA_SMFLUX  && \
      !defined BULK_FLUXES && !defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idUsms, ncFRCid(idUsms,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % sustrG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVsms, ncFRCid(idVsms,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % svstrG(LBi,LBj,1))
#  endif

#  ifdef NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idUsms, ncBLKid(ng), 1,                &
     &                 BLKname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % sustrG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVsms, ncBLKid(ng), 1,                &
     &                 BLKname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % svstrG(LBi,LBj,1))
#  endif
# endif

# if !defined ANA_PAIR && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM       || defined ATM_PRESS)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idPair, ncFRCid(idPair,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 FORCES(ng) % PairG(LBi,LBj,1))
# endif

# if !defined ANA_WWAVE && defined WAVE_DATA
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
#  ifdef WAVES_DIR
      CALL get_2dfldr (ng, iADM, idWdir, ncFRCid(idWdir,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % DwaveG(LBi,LBj,1))
#  endif
#  ifdef WAVES_HEIGHT
      CALL get_2dfldr (ng, iADM, idWamp, ncFRCid(idWamp,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % HwaveG(LBi,LBj,1))
#  endif
#  ifdef WAVES_LENGTH
      CALL get_2dfldr (ng, iADM, idWlen, ncFRCid(idWlen,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % LwaveG(LBi,LBj,1))
#  endif
#  ifdef WAVES_TOP_PERIOD
      CALL get_2dfldr (ng, iADM, idWptp, ncFRCid(idWptp,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % Pwave_topG(LBi,LBj,1))
#  endif
#  ifdef WAVES_BOT_PERIOD
      CALL get_2dfldr (ng, iADM, idWpbt, ncFRCid(idWpbt,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % Pwave_botG(LBi,LBj,1))
#  endif
#  if defined WAVES_UB
      CALL get_2dfldr (ng, iADM, idWorb, ncFRCid(idWorb,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % Ub_swanG(LBi,LBj,1))
#  endif
#  if defined TKE_WAVEDISS
      CALL get_2dfldr (ng, iADM, idWdis, ncFRCid(idWdis,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % Wave_dissipG(LBi,LBj,1))
#  endif
#  if defined SVENDSEN_ROLLER
      CALL get_2dfldr (ng, iADM, idWbrk, ncFRCid(idWbrk,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % Wave_breakG(LBi,LBj,1))
#  endif
# endif

# ifdef SOLVE3D

#  if !defined ANA_CLOUD && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idCfra, ncFRCid(idCfra,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % cloudG(LBi,LBj,1))
#  endif

#  if !defined ANA_SRFLUX && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idSrad, ncFRCid(idSrad,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % srflxG(LBi,LBj,1))
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      !defined LONGWAVE    && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idLrad, ncFRCid(idLrad,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % lrflxG(LBi,LBj,1))
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idLdwn, ncFRCid(idLdwn,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % lrflxG(LBi,LBj,1))
#  endif

#  if !defined ANA_TAIR && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM       || \
      (defined SHORTWAVE   && defined ANA_SRFLUX && defined ALBEDO))
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idTair, ncFRCid(idTair,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % TairG(LBi,LBj,1))
#  endif

#  if !defined ANA_HUMIDITY && \
      ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
       defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idQair, ncFRCid(idQair,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % HairG(LBi,LBj,1))
#  endif

#  if !defined ANA_RAIN    && \
      (defined BULK_FLUXES && !defined NL_BULK_FLUXES)
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idrain, ncFRCid(idrain,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % rainG(LBi,LBj,1))
#  endif

#  if !defined ANA_STFLUX && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idTsur(itemp),                         &
     &                 ncFRCid(idTsur(itemp),ng),                       &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % stflxG(LBi,LBj,1,itemp))
#  endif

#  ifdef NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idTsur(itemp), ncBLKid(ng), 1,         &
     &                 BLKname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % stflxG(LBi,LBj,1,itemp))
#  endif

#  if !defined ANA_SST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST) and heat flux sensitivity to SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idSSTc, ncFRCid(idSSTc,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                FORCES(ng) % sstG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, iddQdT, ncFRCid(iddQdT,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % dqdtG(LBi,LBj,1))
#  endif

#  ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idTbot(itemp),                         &
     &                 ncFRCid(idTbot(itemp),ng),                       &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 FORCES(ng) % btflxG(LBi,LBj,1,itemp))
#  endif

#  ifdef SALINITY
#   if !(defined ANA_SSFLUX || defined EMINUSP || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idsfwf, ncFRCid(idsfwf,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 FORCES(ng) % stflxG(LBi,LBj,1,isalt))
#   endif

#   if defined EMINUSP && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idEmPf, ncBLKid(ng), 1,                &
     &                 BLKname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 FORCES(ng) % stflxG(LBi,LBj,1,isalt))
#   endif

#   if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idSSSc, ncFRCid(idSSSc,ng),            &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 FORCES(ng) % sssG(LBi,LBj,1))
#   endif

#   ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfldr (ng, iADM, idTbot(isalt),                         &
     &                 ncFRCid(idTbot(isalt),ng),                       &
     &                 nFfiles(ng), FRCname(1,ng), update(1),           &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 FORCES(ng) % btflxG(LBi,LBj,1,isalt))
#   endif
#  endif

#  if defined SEDIMENT_NOT_YET || defined BIOLOGY
#   ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfldr (ng, iADM, idTsur(i), ncFRCid(idTsur(i),ng),    &
     &                   nFfiles(ng), FRCname(1,ng), update(1),         &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   FORCES(ng) % stflxG(LBi,LBj,1,i))
      END DO
#   endif

#   ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfldr (ng, iADM, idTbot(i), ncFRCid(idTbot(i),ng),    &
     &                   nFfiles(ng), FRCname(1,ng), update(1),         &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   FORCES(ng) % btflxG(LBi,LBj,1,i))
      END DO
#   endif
#  endif
# endif

# ifdef OBC_DATA
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================
!
#  ifndef ANA_FSOBC
#   ifdef WEST_FSOBC
      CALL get_ngfldr (ng, iADM, idZbry(iwest), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % zetaG_west(JLB,1))
#   endif
#   ifdef EAST_FSOBC
      CALL get_ngfldr (ng, iADM, idZbry(ieast), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % zetaG_east(JLB,1))
#   endif
#   ifdef SOUTH_FSOBC
      CALL get_ngfldr (ng, iADM, idZbry(isouth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % zetaG_south(ILB,1))
#   endif
#   ifdef NORTH_FSOBC
      CALL get_ngfldr (ng, iADM, idZbry(inorth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % zetaG_north(ILB,1))
#   endif
#  endif
#  ifndef ANA_M2OBC
#   ifdef WEST_M2OBC
      CALL get_ngfldr (ng, iADM, idU2bc(iwest), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % ubarG_west(JLB,1))
      CALL get_ngfldr (ng, iADM, idV2bc(iwest), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % vbarG_west(JLB,1))
#   endif
#   ifdef EAST_M2OBC
      CALL get_ngfldr (ng, iADM, idU2bc(ieast), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % ubarG_east(JLB,1))
      CALL get_ngfldr (ng, iADM, idV2bc(ieast), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % vbarG_east(JLB,1))
#   endif
#   ifdef SOUTH_M2OBC
      CALL get_ngfldr (ng, iADM, idU2bc(isouth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % ubarG_south(ILB,1))
      CALL get_ngfldr (ng, iADM, idV2bc(isouth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % vbarG_south(ILB,1))
#   endif
#   ifdef NORTH_M2OBC
      CALL get_ngfldr (ng, iADM, idU2bc(inorth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % ubarG_north(ILB,1))
      CALL get_ngfldr (ng, iADM, idV2bc(inorth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                  &
     &                 BOUNDARY(ng) % vbarG_north(ILB,1))
#   endif
#  endif
#  ifdef SOLVE3D
#   ifndef ANA_M3OBC
#    ifdef WEST_M3OBC
      CALL get_ngfldr (ng, iADM, idU3bc(iwest), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % uG_west(JLB,1,1))
      CALL get_ngfldr (ng, iADM, idV3bc(iwest), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % vG_west(JLB,1,1))
#    endif
#    ifdef EAST_M3OBC
      CALL get_ngfldr (ng, iADM, idU3bc(ieast), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % uG_east(JLB,1,1))
      CALL get_ngfldr (ng, iADM, idV3bc(ieast), ncBRYid(ng), 1,         &
     &                 BRYname(ng), update(1),                          &
     &                 JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % vG_east(JLB,1,1))
#    endif
#    ifdef SOUTH_M3OBC
      CALL get_ngfldr (ng, iADM, idU3bc(isouth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % uG_south(ILB,1,1))
      CALL get_ngfldr (ng, iADM, idV3bc(isouth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % vG_south(ILB,1,1))
#    endif
#    ifdef NORTH_M3OBC
      CALL get_ngfldr (ng, iADM, idU3bc(inorth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % uG_north(ILB,1,1))
      CALL get_ngfldr (ng, iADM, idV3bc(inorth), ncBRYid(ng), 1,        &
     &                 BRYname(ng), update(1),                          &
     &                 ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),          &
     &                 BOUNDARY(ng) % vG_north(ILB,1,1))
#    endif
#   endif
#   ifndef ANA_TOBC
#    ifdef WEST_TOBC
      DO i=1,NT(ng)
        CALL get_ngfldr (ng, iADM, idTbry(iwest,i), ncBRYid(ng), 1,     &
     &                   BRYname(ng), update(1),                        &
     &                   JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),        &
     &                   BOUNDARY(ng) % tG_west(JLB,1,1,i))
      END DO
#    endif
#    ifdef EAST_TOBC
      DO i=1,NT(ng)
        CALL get_ngfldr (ng, iADM, idTbry(ieast,i), ncBRYid(ng), 1,     &
     &                   BRYname(ng), update(1),                        &
     &                   JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),        &
     &                   BOUNDARY(ng) % tG_east(JLB,1,1,i))
      END DO
#    endif
#    ifdef SOUTH_TOBC
      DO i=1,NT(ng)
        CALL get_ngfldr (ng, iADM, idTbry(isouth,i), ncBRYid(ng), 1,    &
     &                   BRYname(ng), update(1),                        &
     &                   ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),        &
     &                   BOUNDARY(ng) % tG_south(ILB,1,1,i))
      END DO
#    endif
#    ifdef NORTH_TOBC
      DO i=1,NT(ng)
        CALL get_ngfldr (ng, iADM, idTbry(inorth,i), ncBRYid(ng), 1,    &
     &                   BRYname(ng), update(1),                        &
     &                   ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),        &
     &                   BOUNDARY(ng) % tG_north(ILB,1,1,i))
      END DO
#    endif
#   endif
#  endif
# endif

# ifdef CLM_FILE
!
!=======================================================================
!  Read in climatology data from CLIMATOLOGY NetCDF file.
!=======================================================================
!
#  if !defined ANA_SSH && defined ZCLIMATOLOGY
      CALL get_2dfldr (ng, iADM, idSSHc, ncCLMid(ng), 1,                &
     &                 CLMname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 CLIMA(ng) % sshG(LBi,LBj,1))
#  endif
#  if !defined ANA_M2CLIMA && defined M2CLIMATOLOGY
      CALL get_2dfldr (ng, iADM, idUbcl, ncCLMid(ng), 1,                &
     &                 CLMname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 CLIMA(ng) % ubarclmG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVbcl, ncCLMid(ng), 1,                &
     &                 CLMname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 CLIMA(ng) % vbarclmG(LBi,LBj,1))
#  endif
#  ifdef SOLVE3D
#   if !defined ANA_TCLIMA && defined TCLIMATOLOGY
      DO i=1,NAT
        CALL get_3dfldr (ng, iADM, idTclm(i), ncCLMid(ng), 1,           &
     &                   CLMname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   CLIMA(ng) % tclmG(LBi,LBj,1,1,i))
      END DO
#   endif
#   if !defined ANA_M3CLIMA && defined M3CLIMATOLOGY
      CALL get_3dfldr (ng, iADM, idUclm, ncCLMid(ng), 1,                &
     &                 CLMname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#    endif
     &                 CLIMA(ng) % uclmG(LBi,LBj,1,1))
      CALL get_3dfldr (ng, iADM, idVclm, ncCLMid(ng), 1,                &
     &                 CLMname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#    endif
     &                 CLIMA(ng) % vclmG(LBi,LBj,1,1))
#   endif
#  endif
# endif

# if defined ASSIMILATION || defined NUDGING
!
!=======================================================================
!  Read in assimilation data from ASSIMILATION NetCDF files.
!=======================================================================
!
#  ifdef ASSIMILATION_SSH
!
!-----------------------------------------------------------------------
!  Read in sea surface height observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SSH(ng)) THEN
        CALL get_2dfldr (ng, iADM, idSSHo, ncSSHid(ng), 1,              &
     &                   SSHname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   OBS(ng) % SSHobs(LBi,LBj))
        CALL get_2dfldr (ng, iADM, idSSHe, ncSSHid(ng), 1,              &
     &                   SSHname(ng), update(2),                        &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   OBS(ng) % EobsSSH(LBi,LBj))
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
#  endif

#  if defined NUDGING_SST || defined ASSIMILATION_SST
!
!-----------------------------------------------------------------------
!  Read in sea surface temperature observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SST(ng)) THEN
#   ifdef NUDGING_SST
        CALL get_2dfldr (ng, iADM, idSSTo, ncSSTid(ng), 1,              &
     &                   SSTname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % SSTdat(LBi,LBj,1))
        CALL get_2dfldr (ng, iADM, idSSTe, ncSSTid(ng), 1,              &
     &                   SSTname(ng), update(2),                        &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % EdatSST(LBi,LBj,1))
#   else
        CALL get_2dfldr (ng, iADM, idSSTo, ncSSTid(ng), 1,              &
     &                   SSTname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % SSTobs(LBi,LBj))
        CALL get_2dfldr (ng, iADM, idSSTe, ncSSTid(ng), 1,              &
     &                   SSTname(ng), update(2),                        &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % EobsSST(LBi,LBj))
#   endif
        IF (update(1).and.update(2)) update_SST(ng)=.TRUE.
      END IF
#  endif

#  if defined NUDGING_T || defined ASSIMILATION_T
!
!-----------------------------------------------------------------------
!  Read in tracers observations and error variance.
!-----------------------------------------------------------------------
!
      DO i=1,NAT
        IF (assi_T(i,ng)) THEN
#   ifdef NUDGING_T
          CALL get_3dfldr (ng, iADM, idTobs(i), ncTOBSid(ng), 1,        &
     &                     TOBSname(ng), update(1),                     &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,          &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     OBS(ng) % Tdat(LBi,LBj,1,1,i))
          CALL get_3dfldr (ng, iADM, idTerr(i), ncTOBSid(ng), 1,        &
     &                     TOBSname(ng), update(2),                     &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,          &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     OBS(ng) % EdatT(LBi,LBj,1,1,i))
          IF (update(1).and.update(2)) update_T(i,ng)=.TRUE.
#   else
          CALL get_3dfldr (ng, iADM, idTobs(i), ncTOBSid(ng), 1,        &
     &                     TOBSname(ng), update(1),                     &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,          &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     OBS(ng) % Tobs(LBi,LBj,1,i))
          CALL get_3dfldr (ng, iADM, idTerr(i), ncTOBSid(ng), 1,        &
     &                     TOBSname(ng), update(2),                     &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,          &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     OBS(ng) % EobsT(LBi,LBj,1,i))
          IF (update(1).and.update(2)) THEN
            update_T(i,ng)=.TRUE.
            tTobs(1,i,ng)=Vtime(1,idTobs(i),ng)
            tsTobs(i,ng)=Vtime(1,idTobs(i),ng)*day2sec
            Finfo(7,idTobs(i),ng)=tsTobs(i,ng)
            Finfo(7,idTerr(i),ng)=tsTobs(i,ng)
            EobsTmin(i,ng)=Finfo(8,idTerr(i),ng)
            EobsTmax(i,ng)=Finfo(9,idTerr(i),ng)
          END IF
#   endif
        END IF
      END DO
#  endif

#  if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
!
!-----------------------------------------------------------------------
!  Read in surface current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UVsur(ng)) THEN
#   ifdef NUDGING_UVsur
        CALL get_2dfldr (ng, iADM, idUsur, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(1),                       &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Usurdat(LBi,LBj,1))
        CALL get_2dfldr (ng, iADM, idVsur, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(2),                       &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Vsurdat(LBi,LBj,1))
        CALL get_2dfldr (ng, iADM, idUVse, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(3),                       &
     &                   LBi, UBi, LBj, UBj, 2, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % EdatVsur(LBi,LBj,1))
#   else
        CALL get_2dfldr (ng, iADM, idUsur, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(1),                       &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Usur(LBi,LBj))
        CALL get_2dfldr (ng, iADM, idVsur, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(2),                       &
     &                   LBi, UBi, LBj, UBj, 1 ,1,                      &
#    ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Vsur(LBi,LBj))
        CALL get_2dfldr (ng, iADM, idUVse, ncVSURid(ng), 1,             &
     &                   VSURname(ng), update(3),                       &
     &                   LBi, UBi, LBj, UBj, 1, 1,                      &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                   OBS(ng) % EobsVsur(LBi,LBj))
#   endif
        IF (update(1).and.update(2).and.update(3)) update_UVsur=.TRUE.
      END IF
#  endif

#  if defined NUDGING_UV || defined ASSIMILATION_UV
!
!-----------------------------------------------------------------------
!  Read in horizontal current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UV(ng)) THEN
#   ifdef NUDGING_UV
        CALL get_3dfldr (ng, iADM, idUobs, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(1),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Udat(LBi,LBj,1,1))
        CALL get_3dfldr (ng, iADM, idVobs, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(2),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Vdat(LBi,LBj,1,1))
        CALL get_3dfldr (ng, iADM, idUVer, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(3),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % EdatUV(LBi,LBj,1,1))
        IF (update(1).and.update(2).and.update(3)) update_UV(ng)=.TRUE.
#   else
        CALL get_3dfldr (ng, iADM, idUobs, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(1),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Uobs(LBi,LBj,1))
        CALL get_3dfldr (ng, iADM, idVobs, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(2),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % Vobs(LBi,LBj,1))
        CALL get_3dfldr (ng, iADM, idUVer, ncVOBSid(ng), 1,             &
     &                   VOBSname(ng), update(3),                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   OBS(ng) % EobsUV(LBi,LBj,1))
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
#   endif
      END IF
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
      CALL get_2dfldr (ng, iADM, idFsur, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 OCEAN(ng) % zetaG(LBi,LBj,1))
!
!  Read 2D momentum.
!
      CALL get_2dfldr (ng, iADM, idUbar, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#  endif
     &                 OCEAN(ng) % ubarG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVbar, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#  ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#  endif
     &                 OCEAN(ng) % vbarG(LBi,LBj,1))
#  ifdef FORWARD_RHS
!
!  Read in variables associated with 2D right-hand-side terms.
!
      CALL get_2dfldr (ng, iADM, idRzet, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 OCEAN(ng) % rzetaG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idRu2d, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 OCEAN(ng) % rubarG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idRv2d, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 OCEAN(ng) % rvbarG(LBi,LBj,1))
#  endif
#  ifdef SOLVE3D
!
!  Read in variables associated with time-averaged 2D momentum terms.
!
      CALL get_2dfldr (ng, iADM, idUfx1, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 COUPLING(ng) % DU_avg1G(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idUfx2, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 COUPLING(ng) % DU_avg2G(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVfx1, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 COUPLING(ng) % DV_avg1G(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idVfx2, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 COUPLING(ng) % DV_avg2G(LBi,LBj,1))
!
!  Read in 3D momentum.
!
      CALL get_3dfldr (ng, iADM, idUvel, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,              &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 OCEAN(ng) % uG(LBi,LBj,1,1))
      CALL get_3dfldr (ng, iADM, idVvel, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,              &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 OCEAN(ng) % vG(LBi,LBj,1,1))
#   ifdef FORWARD_RHS
!
!  Read in variables associated with 3D momentum right-hand-side terms.
!
      CALL get_2dfldr (ng, iADM, idRuct, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#    endif
     &                 COUPLING(ng) % rufrcG(LBi,LBj,1))
      CALL get_2dfldr (ng, iADM, idRvct, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#    endif
     &                 COUPLING(ng) % rvfrcG(LBi,LBj,1))
      CALL get_3dfldr (ng, iADM, idRu3d, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng),  2, 1,             &
#    ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#    endif
     &                 OCEAN(ng) % ruG(LBi,LBj,1,1))
      CALL get_3dfldr (ng, iADM, idRv3d, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#    endif
     &                 OCEAN(ng) % rvG(LBi,LBj,1,1))
#   endif
!
!  Read in 3D tracers.
!
      DO i=1,NT(ng)
        CALL get_3dfldr (ng, iADM, idTvar(i), ncFWDid(ng), 1,           &
     &                   FWDname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,            &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   OCEAN(ng) % tG(LBi,LBj,1,1,i))
      END DO

#   ifdef FORWARD_MIXING
!
! Read in vertical mixing variables.
!
      DO i=1,NAT
        scale=Fscale(idDiff(i),ng)                    ! save and rescale
        Fscale(idDiff(i),ng)=ad_Akt_fac(i,ng)
        CALL get_3dfldr (ng, iADM, idDiff(i), ncFWDid(ng), 1,           &
     &                   FWDname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   MIXING(ng) % AktG(LBi,LBj,0,1,i))
        Fscale(idDiff(i),ng)=scale
      END DO

      scale=Fscale(idVvis,ng)                         ! save and rescale
      Fscale(idVvis,ng)=ad_Akv_fac(ng)
      CALL get_3dfldr (ng, iADM, idVvis, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % AkvG(LBi,LBj,0,1))
      Fscale(idVvis,ng)=scale
#   endif

#   if defined MY25_MIXING_NOT_YET || defined GLS_MIXING_NOT_YET
!
!  Read in turbulent kinetic energy.
!
      CALL get_3dfldr (ng, iADM, idMtke, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % tkeG(LBi,LBj,0,1))
!
!  Read in turbulent kinetic energy times length scale.
!
      CALL get_3dfldr (ng, iADM, idMtls, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % glsG(LBi,LBj,0,1))
!
!  Read in vertical mixing length scale.
!
      CALL get_3dfldr (ng, iADM, idVmLS, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % LscaleG(LBi,LBj,0,1))
!
!  Read in vertical mixing coefficient for turbulent kinetic energy.
!
      CALL get_3dfldr (ng, iADM, idVmKK, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % AkkG(LBi,LBj,0,1))
#    ifdef GLS_MIXING_NOT_YET
!
!  Read in vertical mixing coefficient for turbulent length scale.
!
      CALL get_3dfldr (ng, iADM, idVmKP, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,              &
#     ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#     endif
     &                 MIXING(ng) % AkpG(LBi,LBj,0,1))
#    endif
#   endif

#   ifdef LMD_MIXING_NOT_YET
!
!  Read in depth of surface oceanic boundary layer.
!
      CALL get_2dfldr (ng, iADM, idHsbl, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % hsblG(LBi,LBj,1))
#   endif
#   ifdef LMD_BKPP_NOT_YET
!
!  Read in depth of bottom oceanic boundary layer.
!
      CALL get_2dfldr (ng, iADM, idHbbl, ncFWDid(ng), 1,                &
     &                 FWDname(ng), update(1),                          &
     &                 LBi, UBi, LBj, UBj, 2, 1,                        &
#    ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#    endif
     &                 MIXING(ng) % hbblG(LBi,LBj,1))
#   endif
#   ifdef LMD_NONLOCAL_NOT_YET
!
!  Read in boundary layer nonlocal transport.
!
      DO i=1,NAT
        CALL get_3dfldr (ng, iADM, idGhat(i), ncFWDid(ng), 1,           &
     &                   FWDname(ng), update(1),                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,            &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   MIXING(ng) % ghatsG(LBi,LBj,0,i))
      END DO
#   endif
#  endif
# endif

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Read in adjoint sensitivity functional.
!-----------------------------------------------------------------------
!
!  Read in free-surface.
!
#  ifdef SENSITIVITY_4DVAR
      IF (LsenPSAS(ng)) THEN
#  endif
        IF (SCALARS(ng)%Lstate(isFsur)) THEN
          CALL get_2dfldr (ng, iADM, idZads, ncADSid(ng), 1,            &
     &                     ADSname(ng), update(1),                      &
     &                     LBi, UBi, LBj, UBj, 2, 1,                    &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     CLIMA(ng) % zeta_adsG(LBi,LBj,1))
        END IF
!
!  Read 2D momentum.
!
        IF (SCALARS(ng)%Lstate(isUbar)) THEN
          CALL get_2dfldr (ng, iADM, idUbas, ncADSid(ng), 1,            &
     &                     ADSname(ng), update(1),                      &
     &                     LBi, UBi, LBj, UBj, 2, 1,                    &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     CLIMA(ng) % ubar_adsG(LBi,LBj,1))
        END IF
!
        IF (SCALARS(ng)%Lstate(isVbar)) THEN
          CALL get_2dfldr (ng, iADM, idVbas, ncADSid(ng), 1,            &
     &                     ADSname(ng), update(1),                      &
     &                     LBi, UBi, LBj, UBj, 2, 1,                    &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     CLIMA(ng) % vbar_adsG(LBi,LBj,1))
        END IF

#  ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        IF (SCALARS(ng)%Lstate(isUvel)) THEN
          CALL get_3dfldr (ng, iADM, idUads, ncADSid(ng), 1,            &
     &                     ADSname(ng), update(1),                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,          &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     CLIMA(ng) % u_adsG(LBi,LBj,1,1))
        END IF
!
        IF (SCALARS(ng)%Lstate(isVvel)) THEN
          CALL get_3dfldr (ng, iADM, idVads, ncADSid(ng), 1,            &
     &                     ADSname(ng), update(1),                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,          &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     CLIMA(ng) % v_adsG(LBi,LBj,1,1))
        END IF
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          IF (SCALARS(ng)%Lstate(isTvar(i))) THEN
            CALL get_3dfldr (ng, iADM, idTads(i), ncADSid(ng), 1,       &
     &                       ADSname(ng), update(1),                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#   endif
     &                       CLIMA(ng) % t_adsG(LBi,LBj,1,1,i))
          END IF
        END DO
#  endif
#  ifdef SENSITIVITY_4DVAR
      END IF
#  endif
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iADM, 3)
# endif
      RETURN
      END SUBROUTINE ad_get_data
#else
      SUBROUTINE ad_get_data
      RETURN
      END SUBROUTINE ad_get_data
#endif
