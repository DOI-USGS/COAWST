#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_set_data (ng, tile)
!
!svn $Id: tl_set_data.F 451 2010-03-04 21:43:03Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  assimilation input data. It time-interpolates between snapshots.    !
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
      CALL wclock_on (ng, iTLM, 4)
# endif
      CALL tl_set_data_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 4)
# endif
      RETURN
      END SUBROUTINE tl_set_data
!
!***********************************************************************
      SUBROUTINE tl_set_data_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
# if defined FORWARD_READ && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
# if defined ASSIMILATION || defined NUDGING
      USE mod_obs
# endif
      USE mod_scalars
# if defined TS_PSOURCE || defined UV_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
!
# ifdef ANALYTICAL
      USE analytical_mod
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef SOLVE3D
#  if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
#   if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#   endif
#  endif
# endif
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
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      logical :: update = .FALSE.
# if defined WET_DRY
      logical :: bry_update
# endif

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, itrc, j, k, order

      real(r8) :: Zr, cff, cff1, cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) ::  work1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) ::  work2

# include "set_bounds.h"

# ifdef SOLVE3D

#  ifdef CLOUDS
!
!-----------------------------------------------------------------------
!  Set cloud fraction (nondimensional).  Notice that clouds are
!  processed first in case that they are used to adjust shortwave
!  radiation.
!-----------------------------------------------------------------------
!
#   ifdef ANA_CLOUD
      CALL ana_cloud (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idCfra,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%cloudG,                           &
     &                     FORCES(ng)%cloud,                            &
     &                     update)
#   endif
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM || \
      (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO)
!
!-----------------------------------------------------------------------
!  Set surface air temperature (degC).
!-----------------------------------------------------------------------
!
#   ifdef ANA_TAIR
      CALL ana_tair (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idTair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%TairG,                            &
     &                     FORCES(ng)%Tair,                             &
     &                     update)
#   endif
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM || \
      (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO)
!
!-----------------------------------------------------------------------
!  Set surface air relative or specific humidity.
!-----------------------------------------------------------------------
!
#   ifdef ANA_HUMIDITY
      CALL ana_humid (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idQair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HairG,                            &
     &                     FORCES(ng)%Hair,                             &
     &                     update)
#   endif
#  endif

#  ifdef SHORTWAVE
!
!-----------------------------------------------------------------------
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SRFLUX
      CALL ana_srflux (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
#   endif
#   ifdef DIURNAL_SRFLUX
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      CALL ana_srflux (ng, tile, iTLM)
#   endif
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      !defined LONGWAVE    && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iTLM, idLrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
#  endif

#  if defined LONGWAVE_OUT && \
      (defined BULK_FLUXES && !defined NL_BULK_FLUXES)
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iTLM, idLdwn,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM
!
!-----------------------------------------------------------------------
!  Set surface winds (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_WINDS
      CALL ana_winds (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idUair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%UwindG,                           &
     &                     FORCES(ng)%Uwind,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, idVair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%VwindG,                           &
     &                     FORCES(ng)%Vwind,                            &
     &                     update)
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
#     if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%UWind)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%VWind)
#     endif
#     ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      FORCES(ng)%UWind,                           &
     &                      FORCES(ng)%VWind)
#     endif
      END IF
#    endif
#   endif
#  endif

#  if defined BULK_FLUXES && !defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set rain fall rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_RAIN
      CALL ana_rain (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idrain,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%rainG,                            &
     &                     FORCES(ng)%rain,                             &
     &                     update)
#   endif
#  endif

#  if !defined BULK_FLUXES || defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_STFLUX
      CALL ana_stflux (ng, tile, iTLM, itemp)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idTsur(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%stflx (:,:,itemp),                &
     &                     update)
#   endif
#  endif

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Set sea surface temperature (SST) and heat flux sensitivity to
!  SST (dQdSST) which are used for surface heat flux correction.
!-----------------------------------------------------------------------
!
#   ifdef ANA_SST
      CALL ana_sst (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idSSTc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sstG,                             &
     &                     FORCES(ng)%sst,                              &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, iddQdT,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%dqdtG,                            &
     &                     FORCES(ng)%dqdt,                             &
     &                     update)
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_BTFLUX
      CALL ana_btflux (ng, tile, iTLM, itemp)
#  else
      CALL set_2dfld_tile (ng, tile, iTLM, idTbot(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%btflx (:,:,itemp),                &
     &                     update)
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater (E-P) flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SSFLUX
      CALL ana_stflux (ng, tile, iTLM, isalt)
#   else
#    if !(defined EMINUSP || defined SRELAXATION)
      CALL set_2dfld_tile (ng, tile, iTLM, idsfwf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%stflx (:,:,isalt),                &
     &                     update)
#    endif
#    if defined EMINUSP && defined NL_BULK_FLUXES
      CALL set_2dfld_tile (ng, tile, iTLM, idEmPf,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%stflx (:,:,isalt),                &
     &                     update)
#    endif
#   endif

#   if defined SCORRECTION || defined SRELAXATION
!
!-----------------------------------------------------------------------
!  Set surface salinity for freshwater flux correction.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iTLM)
#    else
      CALL set_2dfld_tile (ng, tile, iTLM, idSSSc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssG,                             &
     &                     FORCES(ng)%sss,                              &
     &                     update)
#    endif
#   endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_BSFLUX
      CALL ana_btflux (ng, tile, iTLM, isalt)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idTbot(isalt),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%btflx (:,:,isalt),                &
     &                     update)
#   endif
#  endif

#  if defined SEDIMENT_NOT_YET || defined BIOLOGY
!
!-----------------------------------------------------------------------
!  Set kinematic surface and bottom pasive tracer fluxes (T m/s).
!-----------------------------------------------------------------------
!
      DO itrc=NAT+1,NT(ng)
#   ifdef ANA_SPFLUX
        CALL ana_stflux (ng, tile, iTLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iTLM, idTsur(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%stflx (:,:,itrc),               &
     &                       update)
#   endif
#   ifdef ANA_BPFLUX
        CALL ana_btflux (ng, tile, iTLM, itrc)
#   else
        CALL set_2dfld_tile (ng, tile, iTLM, idTbot(itrc),              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%btflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%btflx (:,:,itrc),               &
     &                       update)
#   endif
      END DO
#  endif
# endif

# ifndef AIR_OCEAN
#  if !defined BULK_FLUXES || defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SMFLUX
      CALL ana_smflux (ng, tile, iTLM)
#   else
      CALL set_2dfld_tile (ng, tile, iTLM, idUsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sustrG,                           &
     &                     FORCES(ng)%sustr,                            &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, idVsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%svstrG,                           &
     &                     FORCES(ng)%svstr,                            &
     &                     update)
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
#     if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%svstr)
#     endif
#     ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      FORCES(ng)%sustr,                           &
     &                      FORCES(ng)%svstr)
#     endif
      END IF
#    endif
#   endif
#  endif
# endif

# if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
     defined ECOSIM       || defined ATM_PRESS
!
!-----------------------------------------------------------------------
!  Set surface air pressure (mb).
!-----------------------------------------------------------------------
!
#  ifdef ANA_PAIR
      CALL ana_pair (ng, tile, iTLM)
#  else
      CALL set_2dfld_tile (ng, tile, iTLM, idPair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PairG,                            &
     &                     FORCES(ng)%Pair,                             &
     &                     update)
#  endif
# endif

# ifdef WAVE_DATA
!
!-----------------------------------------------------------------------
!  Set surface wind-induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
#  ifdef ANA_WWAVE
      CALL ana_wwave (ng, tile, iTLM)
#  else
#   ifdef WAVES_DIR
      CALL set_2dfld_tile (ng, tile, iTLM, idWdir,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%DwaveG,                           &
     &                     FORCES(ng)%Dwave,                            &
     &                     update)
#    ifdef CURVGRID
!
!  If input point-data, rotate direction to curvilinear coordinates.
!
      IF (.not.Linfo(1,idWdir,ng).or.                                   &
     &    (Iinfo(5,idWdir,ng).ne.Lm(ng)+2).or.                          &
     &    (Iinfo(6,idWdir,ng).ne.Mm(ng)+2)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            FORCES(ng)%Dwave(i,j)=FORCES(ng)%Dwave(i,j)-                &
     &                            GRID(ng)%angler(i,j)
          END DO
        END DO
      END IF
#     if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dwave)
#     endif
#     ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Dwave)
#     endif
#    endif
#   endif
#   ifdef WAVES_HEIGHT
      CALL set_2dfld_tile (ng, tile, iTLM, idWamp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HwaveG,                           &
     &                     FORCES(ng)%Hwave,                            &
     &                     update)
#   endif
#   ifdef WAVES_LENGTH
      CALL set_2dfld_tile (ng, tile, iTLM, idWlen,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%LwaveG,                           &
     &                     FORCES(ng)%Lwave,                            &
     &                     update)
#   endif
#   ifdef WAVES_TOP_PERIOD
      CALL set_2dfld_tile (ng, tile, iTLM, idWptp,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_topG,                       &
     &                     FORCES(ng)%Pwave_top,                        &
     &                     update)
#   endif
#   ifdef WAVES_BOT_PERIOD
      CALL set_2dfld_tile (ng, tile, iTLM, idWpbt,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Pwave_botG,                       &
     &                     FORCES(ng)%Pwave_bot,                        &
     &                     update)
#   endif
#   if defined WAVES_UB
      CALL set_2dfld_tile (ng, tile, iTLM, idWorb,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Ub_swanG,                         &
     &                     FORCES(ng)%Ub_swan,                          &
     &                     update)
#   endif
#   if defined TKE_WAVEDISS
      CALL set_2dfld_tile (ng, tile, iTLM, idWdis,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Wave_dissipG,                     &
     &                     FORCES(ng)%Wave_dissip,                      &
     &                     update)
#   endif
#   if defined SVENDSEN_ROLLER
      CALL set_2dfld_tile (ng, tile, iTLM, idWbrk,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Wave_breakG,                      &
     &                     FORCES(ng)%Wave_break,                       &
     &                     update)
#   endif
#  endif
# endif

# if defined ECOSIM && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute spectral irradiance and cosine of average zenith angle of
!  downwelling spectral photons.
!-----------------------------------------------------------------------
!
      CALL ana_specir (ng, tile, iTLM)
# endif

# ifdef ANA_SPINNING
!
!-----------------------------------------------------------------------
!  Set time-varying rotation force (centripetal accelerations) for
!  polar coordinate grids.
!-----------------------------------------------------------------------
!
      CALL ana_spinning (ng, tile, iTLM)
# endif

# if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
!
#  ifdef ANA_PSOURCE
        CALL ana_psource (ng, tile iTLM)
#  else
      IF (SOUTH_WEST_TEST) THEN
#   if defined UV_PSOURCE || defined Q_PSOURCE
        CALL set_ngfld (ng, iTLM, idRtra, 1, Nsrc(ng), 1,               &
     &                  1, Nsrc(ng), 1,                                 &
     &                  SOURCES(ng) % QbarG,                            &
     &                  SOURCES(ng) % Qbar,                             &
     &                  update)
#    ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=1,Nsrc(ng)
            SOURCES(ng)%Qsrc(i,k)=SOURCES(ng)%Qbar(i)*                  &
     &                            SOURCES(ng)%Qshape(i,k)
          END DO
        END DO
#    endif
#   endif
#   if defined TS_PSOURCE && defined SOLVE3D
        DO itrc=1,NT(ng)
          IF (LtracerSrc(itrc,ng)) THEN
            CALL set_ngfld (ng, iTLM, idRtrc(itrc), 1, Nsrc(ng), N(ng), &
     &                      1, Nsrc(ng), N(ng),                         &
     &                      SOURCES(ng) % TsrcG(1,1,1,itrc),            &
     &                      SOURCES(ng) % Tsrc(1,1,itrc),               &
     &                      update)
          END IF
        END DO
#   endif
      END IF
#  endif
# endif

# ifdef OBC
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.
!-----------------------------------------------------------------------
!
!  Lower and upper bounds for nontiled arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J

#  ifdef ANA_FSOBC
      CALL ana_fsobc (ng, tile, iTLM)
#  else
      IF (SOUTH_WEST_TEST) THEN
#   ifdef WEST_FSOBC
        CALL set_ngfld (ng, iTLM, idZbry(iwest), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_west(JLB,1),               &
     &                  BOUNDARY(ng) % zeta_west(JLB),                  &
     &                  update)
#   endif
#   ifdef EAST_FSOBC
        CALL set_ngfld (ng, iTLM, idZbry(ieast), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_east(JLB,1),               &
     &                  BOUNDARY(ng) % zeta_east(JLB),                  &
     &                  update)
#   endif
#   ifdef SOUTH_FSOBC
        CALL set_ngfld (ng, iTLM, idZbry(isouth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1 ,1,                                 &
     &                  BOUNDARY(ng) % zetaG_south(ILB,1),              &
     &                  BOUNDARY(ng) % zeta_south(ILB),                 &
     &                  update)
#   endif
#   ifdef NORTH_FSOBC
        CALL set_ngfld (ng, iTLM, idZbry(inorth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_north(ILB,1),              &
     &                  BOUNDARY(ng) % zeta_north(ILB),                 &
     &                  update)
#   endif
      END IF
#  endif

#  if defined WET_DRY
!
! Ensure that water level on boundary cells is above bed elevation.
!
#   ifdef WEST_FSOBC
      bry_update=.FALSE.
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=Dcrit(ng)-GRID(ng)%h(0,j)
          IF (BOUNDARY(ng)%zeta_west(j).le.cff) THEN
            BOUNDARY(ng)%zeta_west(j)=cff
          END IF
        END DO
        bry_update=.TRUE.
      END IF
#    ifdef DISTRIBUTE
      CALL mp_boundary (ng, iTLM, JstrR, JendR, JLB, JUB, 1, 1,         &
     &                  bry_update,                                     &
     &                  BOUNDARY(ng)%zeta_west)
#    endif
#   endif

#   ifdef EAST_FSOBC
      bry_update=.FALSE.
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=Dcrit(ng)-GRID(ng)%h(Lm(ng)+1,j)
          IF (BOUNDARY(ng)%zeta_east(j).le.cff) THEN
            BOUNDARY(ng)%zeta_east(j)=cff
          END IF
        END DO
        bry_update=.TRUE.
      END IF
#    ifdef DISTRIBUTE
      CALL mp_boundary (ng, iTLM, JstrR, JendR, JLB, JUB, 1, 1,         &
     &                  bry_update,                                     &
     &                  BOUNDARY(ng)%zeta_east)
#    endif
#   endif

#   ifdef SOUTH_FSOBC
      bry_update=.FALSE.
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          cff=Dcrit(ng)-GRID(ng)%h(i,0)
          IF (BOUNDARY(ng)%zeta_south(i).le.cff) THEN
            BOUNDARY(ng)%zeta_south(i)=cff
          END IF
        END DO
        bry_update=.TRUE.
      END IF
#    ifdef DISTRIBUTE
      CALL mp_boundary (ng, iTLM, IstrR, IendR, ILB, IUB, 1, 1,         &
     &                  bry_update,                                     &
     &                  BOUNDARY(ng)%zeta_south)
#    endif
#   endif

#   ifdef NORTH_FSOBC
      bry_update=.FALSE.
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          cff=Dcrit(ng)-GRID(ng)%h(i,Mm(ng)+1)
          IF (BOUNDARY(ng)%zeta_north(i).le.cff) THEN
            BOUNDARY(ng)%zeta_north(i)=cff
          END IF
        END DO
        bry_update=.TRUE.
      END IF
#    ifdef DISTRIBUTE
      CALL mp_boundary (ng, iTLM, IstrR, IendR, ILB, IUB, 1, 1,         &
     &                  bry_update,                                     &
     &                  BOUNDARY(ng)%zeta_north)
#    endif
#   endif
#  endif

#  ifdef ANA_M2OBC
      CALL ana_m2obc (ng, tile, iTLM)
#  else
      IF (SOUTH_WEST_TEST) THEN
#   ifdef WEST_M2OBC
        CALL set_ngfld (ng, iTLM, idU2bc(iwest), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_west(JLB,1),               &
     &                  BOUNDARY(ng) % ubar_west(JLB),                  &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV2bc(iwest), JLB, JUB, 1,           &
     &                  1, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_west(JLB,1),               &
     &                  BOUNDARY(ng) % vbar_west(JLB),                  &
     &                  update)
#   endif
#   ifdef EAST_M2OBC
        CALL set_ngfld (ng, iTLM, idU2bc(ieast), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_east(JLB,1),               &
     &                  BOUNDARY(ng) % ubar_east(JLB),                  &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV2bc(ieast), JLB, JUB, 1,           &
     &                  1, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_east(JLB,1),               &
     &                  BOUNDARY(ng) % vbar_east(JLB),                  &
     &                  update)
#   endif
#   ifdef SOUTH_M2OBC
        CALL set_ngfld (ng, iTLM, idU2bc(isouth), ILB, IUB, 1,          &
     &                  1, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_south(ILB,1),              &
     &                  BOUNDARY(ng) % ubar_south(ILB),                 &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV2bc(isouth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_south(ILB,1),              &
     &                  BOUNDARY(ng) % vbar_south(ILB),                 &
     &                  update)
#   endif
#   ifdef NORTH_M2OBC
        CALL set_ngfld (ng, iTLM, idU2bc(inorth), ILB, IUB, 1,          &
     &                  1, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_north(ILB,1),              &
     &                  BOUNDARY(ng) % ubar_north(ILB),                 &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV2bc(inorth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_north(ILB,1),              &
     &                  BOUNDARY(ng) % vbar_north(ILB),                 &
     &                  update)
#   endif
      END IF
#  endif
#  ifdef SOLVE3D
#   ifdef ANA_M3OBC
      CALL ana_m3obc (ng, tile, iTLM)
#   else
      IF (SOUTH_WEST_TEST) THEN
#    ifdef WEST_M3OBC
        CALL set_ngfld (ng, iTLM, idU3bc(iwest), JLB, JUB, N(ng),       &
     &                  0, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_west(JLB,1,1),                &
     &                  BOUNDARY(ng) % u_west(JLB,1),                   &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV3bc(iwest), JLB, JUB, N(ng),       &
     &                  1, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_west(JLB,1,1),                &
     &                  BOUNDARY(ng) % v_west(JLB,1),                   &
     &                  update)
#    endif
#    ifdef EAST_M3OBC
        CALL set_ngfld (ng, iTLM, idU3bc(ieast), JLB, JUB, N(ng),       &
     &                  0, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_east(JLB,1,1),                &
     &                  BOUNDARY(ng) % u_east(JLB,1),                   &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV3bc(ieast), JLB, JUB, N(ng),       &
     &                  1, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_east(JLB,1,1),                &
     &                  BOUNDARY(ng) % v_east(JLB,1),                   &
     &                  update)
#    endif
#    ifdef SOUTH_M3OBC
        CALL set_ngfld (ng, iTLM, idU3bc(isouth), ILB, IUB, N(ng),      &
     &                  1, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_south(ILB,1,1),               &
     &                  BOUNDARY(ng) % u_south(ILB,1),                  &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV3bc(isouth), ILB, IUB, N(ng),      &
     &                  0, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_south(ILB,1,1),               &
     &                  BOUNDARY(ng) % v_south(ILB,1),                  &
     &                  update)
#    endif
#    ifdef NORTH_M3OBC
        CALL set_ngfld (ng, iTLM, idU3bc(inorth), ILB, IUB, N(ng),      &
     &                  1, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_north(ILB,1,1),               &
     &                  BOUNDARY(ng) % u_north(ILB,1),                  &
     &                  update)
        CALL set_ngfld (ng, iTLM, idV3bc(inorth), ILB, IUB, N(ng),      &
     &                  0, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_north(ILB,1,1),               &
     &                  BOUNDARY(ng) % v_north(ILB,1),                  &
     &                  update)
#    endif
      END IF
#   endif
#   ifdef ANA_TOBC
      CALL ana_tobc (ng, tile, iTLM)
#   else
      IF (SOUTH_WEST_TEST) THEN
        DO itrc=1,NT(ng)
#    ifdef WEST_TOBC
          CALL set_ngfld (ng, iTLM, idTbry(iwest,itrc),                 &
     &                    JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_west(JLB,1,1,itrc),         &
     &                    BOUNDARY(ng) % t_west(JLB,1,itrc),            &
     &                    update)
#    endif
#    ifdef EAST_TOBC
          CALL set_ngfld (ng, iTLM, idTbry(ieast,itrc),                 &
     &                    JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_east(JLB,1,1,itrc),         &
     &                    BOUNDARY(ng) % t_east(JLB,1,itrc),            &
     &                    update)
#    endif
#    ifdef SOUTH_TOBC
          CALL set_ngfld (ng, iTLM, idTbry(isouth,itrc),                &
     &                    ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_south(ILB,1,1,itrc),        &
     &                    BOUNDARY(ng) % t_south(ILB,1,itrc),           &
     &                    update)
#    endif
#    ifdef NORTH_TOBC
          CALL set_ngfld (ng, iTLM, idTbry(inorth,itrc),                &
     &                    ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_north(ILB,1,1,itrc),        &
     &                    BOUNDARY(ng) % t_north(ILB,1,itrc),           &
     &                    update)
#    endif
        END DO
      END IF
#   endif
#  endif
# endif

# ifdef ZCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set sea surface height climatology (m).
!-----------------------------------------------------------------------
!
#  ifdef ANA_SSH
      CALL ana_ssh (ng, tile, iTLM)
#  else
      CALL set_2dfld_tile (ng, tile, iTLM, idSSHc,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%sshG,                              &
     &                     CLIMA(ng)%ssh,                               &
     &                     update)
#  endif
# endif

# ifdef M2CLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_M2CLIMA
      CALL ana_m2clima (ng, tile, iTLM)
#  else
      CALL set_2dfld_tile (ng, tile, iTLM, idUbcl,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%ubarclmG,                          &
     &                     CLIMA(ng)%ubarclm,                           &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, idVbcl,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%vbarclmG,                          &
     &                     CLIMA(ng)%vbarclm,                           &
     &                     update)
#  endif
# endif

# if defined SOLVE3D && defined TCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#  ifdef ANA_TCLIMA
      CALL ana_tclima (ng, tile, iTLM)
#  else
      DO itrc=1,NAT
        CALL set_3dfld_tile (ng, tile, iTLM, idTclm(itrc),              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%tclmG(:,:,:,:,itrc),             &
     &                       CLIMA(ng)%tclm (:,:,:,itrc),               &
     &                       update)
      END DO
#  endif
# endif

# if defined SOLVE3D && defined M3CLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_M3CLIMA
      CALL ana_m3clima (ng, tile, iTLM)
#  else
      CALL set_3dfld_tile (ng, tile, iTLM, idUclm,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     CLIMA(ng)%uclmG,                             &
     &                     CLIMA(ng)%uclm,                              &
     &                     update)
      CALL set_3dfld_tile (ng, tile, iTLM, idVclm,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     CLIMA(ng)%vclmG,                             &
     &                     CLIMA(ng)%vclm,                              &
     &                     update)
#  endif
# endif

# ifdef SOLVE3D
#  if defined NUDGING_SST || defined ASSIMILATION_SST
!
!-----------------------------------------------------------------------
!  Set sea surface temperature observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SST(ng)) THEN
#   ifdef NUDGING_SST
        CALL set_2dfld_tile (ng, tile, iTLM, idSSTo,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%SSTdat,                            &
     &                       OBS(ng)%SSTobs,                            &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idSSTe,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%EdatSST,                           &
     &                       OBS(ng)%EobsSST,                           &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_SST(ng)=.FALSE.
          update_T(itemp,ng)=.FALSE.
          exit_flag=0
        END IF
#   endif
!
!  Extend sea surface temperature and associated error variance using
!  provided basis function polynomials.
!
        IF (extend_SST(ng).and.update_SST(ng)) THEN
          IF (SOUTH_WEST_TEST) THEN
            update_SST(ng)=.FALSE.
            update_T(itemp,ng)=.TRUE.
#   ifdef ASSIMILATION_SST
            tTobs(1,itemp,ng)=Vtime(1,idSSTo,ng)
            tsTobs(itemp,ng)=Vtime(1,idSSTo,ng)*day2sec
            Finfo(7,idSSTo,ng)=tsTobs(itemp,ng)
            Finfo(7,idSSTe,ng)=tsTobs(itemp,ng)
            EobsTmin(itemp,ng)=Finfo(8,idSSTe,ng)
            EobsTmax(itemp,ng)=Finfo(9,idSSTe,ng)
#   endif
          END IF
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff=perr_SST(npSST(ng),ng)
                cff1=pcoef_SST(npSST(ng),ng)
                DO order=npSST(ng)-1,0,-1
                  cff=Zr*cff+perr_SST(order,ng)
                  cff1=Zr*cff1+pcoef_SST(order,ng)
                END DO
                OBS(ng)%EobsT(i,j,k,itemp)=MIN(1.0_r8,cff*              &
     &                                         OBS(ng)%EobsSST(i,j))
                OBS(ng)%Tobs(i,j,k,itemp)=cff1*OBS(ng)%SSTobs(i,j)
              END DO
            END DO
          END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OBS(ng)%EobsT(:,:,:,itemp))
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OBS(ng)%Tobs(:,:,:,itemp))
#   endif
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iTLM, 2,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        OBS(ng)%EobsT(:,:,:,itemp),               &
     &                        OBS(ng)%Tobs(:,:,:,itemp))
#   endif
        END IF
      END IF
#  endif

#  if defined NUDGING_T
!
!-----------------------------------------------------------------------
!  Set tracers observations and error variance.
!-----------------------------------------------------------------------
!
      DO itrc=1,NAT
        IF (assi_T(itrc,ng)) THEN
          CALL set_3dfld_tile (ng, tile, iTLM, idTobs(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OBS(ng)%Tdat(:,:,:,:,itrc),              &
     &                         OBS(ng)%Tobs(:,:,:,itrc),                &
     &                         update)
          CALL set_3dfld_tile (ng, tile, iTLM, idTerr(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OBS(ng)%EdatT(:,:,:,:,itrc),             &
     &                         OBS(ng)%EobsT(:,:,:,itrc),               &
     &                         update)
          IF (.not.update.and.SOUTH_WEST_TEST) THEN
            update_T(itrc,ng)=.FALSE.
            exit_flag=0
          END IF
        END IF
      END DO
#  endif

#  if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
!
!-----------------------------------------------------------------------
!  Set surface current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UVsur(ng)) THEN
#   ifdef NUDGING_UVsur
        CALL set_2dfld_tile (ng, tile, iTLM, idUsur,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%Usurdat,                           &
     &                       OBS(ng)%Usur,                              &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idVsur,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%Vsurdat,                           &
     &                       OBS(ng)%Vsur,                              &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idUVse,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%EdatVsur,                          &
     &                       OBS(ng)%EobsVsur,                          &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_UVsur(ng)=.FALSE.
          update_UV(ng)=.FALSE.
          exit_flag=0
        END IF
#   endif
!
!  Extend surface currents observations and associated error variance
!  using provided basis function polynomials.
!
        IF (extend_UV(ng).and.update_UVsur(ng)) THEN
          IF (SOUTH_WEST_TEST) THEN
            update_UVsur(ng)=.FALSE.
            update_UV(ng)=.TRUE.
#   ifdef ASSIMILATION_UVsur
            tVobs(1,ng)=Vtime(1,idVsur,ng)
            tsVobs(ng)=Vtime(1,idVsur,ng)*day2sec
            Finfo(7,idUsur,ng)=tsVobs(ng)
            Finfo(7,idVsur,ng)=tsVobs(ng)
            Finfo(7,idUVse,ng)=tsVobs(ng)
            EobsUVmin(ng)=Finfo(8,idUVse,ng)
            EobsUVmax(ng)=Finfo(9,idUVse,ng)
#   endif
          END IF
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff=perr_V(npUV(ng),ng)
                DO order=npUV(ng)-1,0,-1
                  cff=Zr*cff+perr_V(order,ng)
                END DO
                OBS(ng)%EobsUV(i,j,k)=MIN(1.0_r8,cff+                   &
     &                                    OBS(ng)%EobsVsur(i,j))
              END DO
            END DO
            DO j=JstrV-1,Jend
              DO i=IstrU-1,Iend
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff1=pcoef_U(npUV(ng),ng)
                cff2=pcoef_V(npUV(ng),ng)
                DO order=npUV(ng)-1,0,-1
                  cff1=Zr*cff1+pcoef_U(order,ng)
                  cff2=Zr*cff2+pcoef_V(order,ng)
                END DO
                work1(i,j)=cff1*OBS(ng)%Usur(i,j)-                      &
     &                     cff2*OBS(ng)%Vsur(i,j)
                work2(i,j)=cff2*OBS(ng)%Usur(i,j)+                      &
     &                     cff1*OBS(ng)%Vsur(i,j)
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                OBS(ng)%Uobs(i,j,k)=0.5_r8*(work1(i-1,j)+work1(i,j))
              END DO
              IF (j.ge.JstrV) THEN
                DO i=Istr,Iend
                  OBS(ng)%Vobs(i,j,k)=0.5_r8*(work2(i,j-1)+work2(i,j))
                END DO
              END IF
            END DO
          END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OBS(ng)%EobsUV)
#   endif
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        OBS(ng)%EobsUV)
#   endif
        END IF
      END IF
#  endif

#  ifdef NUDGING_UV
!
!-----------------------------------------------------------------------
!  Set horizontal current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UV(ng)) THEN
        CALL set_3dfld_tile (ng, tile, iTLM, idUobs,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%Udat,                              &
     &                       OBS(ng)%Uobs,                              &
     &                       update)
        CALL set_3dfld_tile (ng, tile, iTLM, idVobs,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%Vdat,                              &
     &                       OBS(ng)%Vobs,                              &
     &                       update)
        CALL set_3dfld_tile (ng, tile, iTLM, idUVer,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%EdatUV,                            &
     &                       OBS(ng)%EobsUV,                            &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_UV(ng)=.FALSE.
          exit_flag=0
        END IF
      END IF
#  endif
# endif

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Set forward solution needed by Tangent/Adjoint models.
!-----------------------------------------------------------------------
!
!  Set forward free-surface.
!
      DO k=1,3
        CALL set_2dfld_tile (ng, tile, iTLM, idFsur,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%zetaG,                           &
     &                       OCEAN(ng)%zeta(:,:,k),                     &
     &                       update)
      END DO
#  ifdef SOLVE3D
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          COUPLING(ng)%Zt_avg1(i,j)=OCEAN(ng)%zeta(i,j,1)
        END DO
      END DO
#  endif
!
!  Set forward 2D momentum.
!
      DO k=1,3
        CALL set_2dfld_tile (ng, tile, iTLM, idUbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%ubarG,                           &
     &                       OCEAN(ng)%ubar(:,:,k),                     &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idVbar,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%vbarG,                           &
     &                       OCEAN(ng)%vbar(:,:,k),                     &
     &                       update)
      END DO
#  ifdef FORWARD_RHS
!
!  Set forward variables associated with 2D right-hand-side terms.
!
      DO k=1,2
        CALL set_2dfld_tile (ng, tile, iTLM, idRzet,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%rzetaG,                          &
     &                       OCEAN(ng)%rzeta(:,:,k),                    &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idRu2d,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%rubarG,                          &
     &                       OCEAN(ng)%rubar(:,:,k),                    &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idRv2d,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%rvbarG,                          &
     &                       OCEAN(ng)%rvbar(:,:,k),                    &
     &                       update)
      END DO
#  endif
#  ifdef SOLVE3D
!
!  Set forward time-averaged barotropic fluxes.
!
      CALL set_2dfld_tile (ng, tile, iTLM, idUfx1,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     COUPLING(ng)%DU_avg1G,                       &
     &                     COUPLING(ng)%DU_avg1,                        &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, idUfx2,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     COUPLING(ng)%DU_avg2G,                       &
     &                     COUPLING(ng)%DU_avg2,                        &
     &                     update)

      CALL set_2dfld_tile (ng, tile, iTLM, idVfx1,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     COUPLING(ng)%DV_avg1G,                       &
     &                     COUPLING(ng)%DV_avg1,                        &
     &                     update)
      CALL set_2dfld_tile (ng, tile, iTLM, idVfx2,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     COUPLING(ng)%DV_avg2G,                       &
     &                     COUPLING(ng)%DV_avg2,                        &
     &                     update)
!
!  Set 3D momentum.
!
      DO k=1,2
        CALL set_3dfld_tile (ng, tile, iTLM, idUvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%uG,                              &
     &                       OCEAN(ng)%u(:,:,:,k),                      &
     &                       update)
        CALL set_3dfld_tile (ng, tile, iTLM, idVvel,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%vG,                              &
     &                       OCEAN(ng)%v(:,:,:,k),                      &
     &                       update)
      END DO

#   ifdef FORWARD_RHS
!
!  Set variables associated with 3D momentum right-hand-side terms.
!
      DO k=1,2
        CALL set_3dfld_tile (ng, tile, iTLM, idRu3d,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%ruG,                             &
     &                       OCEAN(ng)%ru(:,:,:,k),                     &
     &                       update)

        CALL set_3dfld_tile (ng, tile, iTLM, idRv3d,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%rvG,                             &
     &                       OCEAN(ng)%rv(:,:,:,k),                     &
     &                       update)
      END DO
#   endif
!
!  Set 3D tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,3
          CALL set_3dfld_tile (ng, tile, iTLM, idTvar(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OCEAN(ng)%tG(:,:,:,:,itrc),              &
     &                         OCEAN(ng)%t(:,:,:,k,itrc),               &
     &                         update)
        END DO
      END DO

#   ifdef FORWARD_MIXING
!
!  Set forward vertical mixing variables.
!
      DO itrc=1,NAT
        CALL set_3dfld_tile (ng, tile, iTLM, idDiff(itrc),              &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       MIXING(ng)%AktG(:,:,:,:,itrc),             &
     &                       MIXING(ng)%Akt(:,:,:,itrc),                &
     &                       update)
      END DO
      CALL set_3dfld_tile (ng, tile, iTLM, idVvis,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     MIXING(ng)%AkvG,                             &
     &                     MIXING(ng)%Akv,                              &
     &                     update)
#   endif

#   if defined MY25_MIXING_NOT_YET || defined GLS_MIXING_NOT_YET
!
!  Set forward turbulent kinetic energy.
!
      DO k=1,3
        CALL set_3dfld_tile (ng, tile, iTLM, idMtke,                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       MIXING(ng)%tkeG,                           &
     &                       MIXING(ng)%tke(:,:,:,k),                   &
     &                       update)
      END DO
!
!  Set forward turbulent kinetic energy times length scale.
!
      DO k=1,3
        CALL set_3dfld_tile (ng, tile, iTLM, idMtls,                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       MIXING(ng)%glsG,                           &
     &                       MIXING(ng)%gls(:,:,:,k),                   &
     &                       update)
      END DO
!
!  Set forward vertical mixing length scale.
!
      CALL set_3dfld_tile (ng, tile, iTLM, idVmLS,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     MIXING(ng)%LscaleG,                          &
     &                     MIXING(ng)%Lscale,                           &
     &                     update)
!
!  Set forward vertical mixing coefficient for turbulent kinetic energy.
!
      CALL set_3dfld_tile (ng, tile, iTLM, idVmKK,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     MIXING(ng)%AkkG,                             &
     &                     MIXING(ng)%Akk,                              &
     &                     update)
#    ifdef GLS_MIXING_NOT_YET
!
!  Set forward vertical mixing coefficient for turbulent length scale.
!
      CALL set_3dfld_tile (ng, tile, iTLM, idVmKP,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     MIXING(ng)%AkpG,                             &
     &                     MIXING(ng)%Akp,                              &
     &                     update)
#    endif
#   endif

#   ifdef LMD_MIXING_NOT_YET
!
!  Set forward depth of surface oceanic boundary layer.
!
      CALL set_2dfld_tile (ng, tile, iTLM, idHsbl,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     MIXING(ng)%hsblG,                            &
     &                     MIXING(ng)%hsbl,                             &
     &                     update)
#   endif
#   ifdef LMD_BKPP_NOT_YET
!
!  Set forward depth of bottom oceanic boundary layer.
!
      CALL set_2dfld_tile (ng, tile, iTLM, idHbbl,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     MIXING(ng)%hbblG,                            &
     &                     MIXING(ng)%hbbl,                             &
     &                     update)
#   endif
#   ifdef LMD_NONLOCAL_NOT_YET
!
!  Set forward boundary layer nonlocal transport.
!
      DO itrc=1,NAT
        CALL set_3dfld_tile (ng, tile, iTLM, idGhat(itrc),              &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       MIXING(ng)%ghatsG(:,:,:,:,itrc),           &
     &                       MIXING(ng)%ghat(:,:,:,itrc),               &
     &                       update)
      END DO
#   endif
#  endif

#  ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Set weak contraint frequent impulse forcing.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse) THEN
!
!  Set free-surface forcing.
!
        CALL set_2dfld_tile (ng, tile, iTLM, idZtlf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%f_zetaG,                         &
     &                       OCEAN(ng)%f_zeta,                          &
     &                       update)

#   ifndef SOLVE3D
!
!  Set 2D momentum forcing.
!
        CALL set_2dfld_tile (ng, tile, iTLM, idUbtf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%f_ubarG,                         &
     &                       OCEAN(ng)%f_ubar,                          &
     &                       update)
        CALL set_2dfld_tile (ng, tile, iTLM, idVbtf,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng)%f_vbarG,                         &
     &                       OCEAN(ng)%f_vbar,                          &
     &                       update)
#   endif

#   ifdef SOLVE3D
!
!  Set 3D momentum forcing.
!
        CALL set_3dfld_tile (ng, tile, iTLM, idUtlf,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%f_uG,                            &
     &                       OCEAN(ng)%f_u,                             &
     &                       update)
        CALL set_3dfld_tile (ng, tile, iTLM, idVtlf,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng)%f_vG,                            &
     &                       OCEAN(ng)%f_v,                             &
     &                       update)
!
!  Set 3D tracers forcing.
!
        DO itrc=1,NT(ng)
          CALL set_3dfld_tile (ng, tile, iTLM, idTtlf(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OCEAN(ng)%f_tG(:,:,:,:,itrc),            &
     &                         OCEAN(ng)%f_t(:,:,:,itrc),               &
     &                         update)
        END DO
#   endif
      END IF
#  endif
# endif

      RETURN
      END SUBROUTINE tl_set_data_tile
#else
      SUBROUTINE tl_set_data
      RETURN
      END SUBROUTINE tl_set_data
#endif
