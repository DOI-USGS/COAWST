#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_set_data (ng, tile)
!
!svn $Id: ad_set_data.F 795 2016-05-11 01:42:43Z arango $
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
      CALL wclock_on (ng, iADM, 4)
# endif
      CALL ad_set_data_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 4)
# endif

      RETURN
      END SUBROUTINE ad_set_data
!
!***********************************************************************
      SUBROUTINE ad_set_data_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
# if defined FORWARD_READ && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
# ifdef SENSITIVITY_4DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
      USE mod_sources
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
      USE set_2dfldr_mod
# ifdef SOLVE3D
      USE set_3dfldr_mod
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
      CALL ana_cloud (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idCfra,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%cloudG,                          &
     &                      FORCES(ng)%cloud,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
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
      CALL ana_tair (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idTair,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%TairG,                           &
     &                      FORCES(ng)%Tair,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
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
      CALL ana_humid (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idQair,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%HairG,                           &
     &                      FORCES(ng)%Hair,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef SHORTWAVE
!
!-----------------------------------------------------------------------
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SRFLUX
      CALL ana_srflux (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idSrad,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%srflxG,                          &
     &                      FORCES(ng)%srflx,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#   ifdef DIURNAL_SRFLUX
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      CALL ana_srflux (ng, tile, iADM)
#   endif
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      !defined LONGWAVE    && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfldr_tile (ng, tile, iADM, idLrad,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%lrflxG,                          &
     &                      FORCES(ng)%lrflx,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined LONGWAVE_OUT && \
      (defined BULK_FLUXES && !defined NL_BULK_FLUXES)
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfldr_tile (ng, tile, iADM, idLdwn,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%lrflxG,                          &
     &                      FORCES(ng)%lrflx,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM
!
!-----------------------------------------------------------------------
!  Set surface winds (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_WINDS
      CALL ana_winds (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idUair,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%UwindG,                          &
     &                      FORCES(ng)%Uwind,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfldr_tile (ng, tile, iADM, idVair,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%VwindG,                          &
     &                      FORCES(ng)%Vwind,                           &
     &                      update)
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
        CALL mp_exchange2d (ng, tile, iADM, 2,                          &
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

#  if defined BULK_FLUXES && !defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set rain fall rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_RAIN
      CALL ana_rain (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idrain,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%rainG,                           &
     &                      FORCES(ng)%rain,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if !defined BULK_FLUXES || defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_STFLUX
      CALL ana_stflux (ng, tile, iADM, itemp)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idTsur(itemp),              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%stflxG(:,:,:,itemp),             &
     &                      FORCES(ng)%stflx (:,:,itemp),               &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
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
      CALL ana_sst (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idSSTc,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%sstG,                            &
     &                      FORCES(ng)%sst,                             &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
!
#   ifdef ANA_DQDSST
      CALL ana_dqdsst (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, iddQdT,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%dqdtG,                           &
     &                      FORCES(ng)%dqdt,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_BTFLUX
      CALL ana_btflux (ng, tile, iADM, itemp)
#  else
      CALL set_2dfldr_tile (ng, tile, iADM, idTbot(itemp),              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%btflxG(:,:,:,itemp),             &
     &                      FORCES(ng)%btflx (:,:,itemp),               &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater (E-P) flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SSFLUX
      CALL ana_stflux (ng, tile, iADM, isalt)
#   else
#    if !(defined EMINUSP || defined SRELAXATION)
      CALL set_2dfldr_tile (ng, tile, iADM, idsfwf,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%stflxG(:,:,:,isalt),             &
     &                      FORCES(ng)%stflx (:,:,isalt),               &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#    endif

#    if defined EMINUSP && defined NL_BULK_FLUXES
      CALL set_2dfldr_tile (ng, tile, iADM, idEmPf,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%stflxG(:,:,:,isalt),             &
     &                      FORCES(ng)%stflx (:,:,isalt),               &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif

#   if defined SCORRECTION || defined SRELAXATION
!
!-----------------------------------------------------------------------
!  Set surface salinity for freshwater flux correction.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iADM)
#    else
      CALL set_2dfldr_tile (ng, tile, iADM, idSSSc,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%sssG,                            &
     &                      FORCES(ng)%sss,                             &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_BSFLUX
      CALL ana_btflux (ng, tile, iADM, isalt)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idTbot(isalt),              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%btflxG(:,:,:,isalt),             &
     &                      FORCES(ng)%btflx (:,:,isalt),               &
     &                      update)
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
        CALL ana_stflux (ng, tile, iADM, itrc)
#   else
        CALL set_2dfldr_tile (ng, tile, iADM, idTsur(itrc),             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%stflxG(:,:,:,itrc),            &
     &                        FORCES(ng)%stflx (:,:,itrc),              &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef ANA_BPFLUX
        CALL ana_btflux (ng, tile, iADM, itrc)
#   else
        CALL set_2dfldr_tile (ng, tile, iADM, idTbot(itrc),             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%btflxG(:,:,:,itrc),            &
     &                        FORCES(ng)%btflx (:,:,itrc),              &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
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
      CALL ana_smflux (ng, tile, iADM)
#   else
      CALL set_2dfldr_tile (ng, tile, iADM, idUsms,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%sustrG,                          &
     &                      FORCES(ng)%sustr,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfldr_tile (ng, tile, iADM, idVsms,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%svstrG,                          &
     &                      FORCES(ng)%svstr,                           &
     &                      update)
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
        CALL mp_exchange2d (ng, tile, iADM, 2,                          &
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

# if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
     defined ECOSIM       || defined ATM_PRESS
!
!-----------------------------------------------------------------------
!  Set surface air pressure (mb).
!-----------------------------------------------------------------------
!
#  ifdef ANA_PAIR
      CALL ana_pair (ng, tile, iADM)
#  else
      SetBC=.TRUE.
!     SetBC=.FALSE.
      CALL set_2dfldr_tile (ng, tile, iADM, idPair,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%PairG,                           &
     &                      FORCES(ng)%Pair,                            &
     &                      update, SetBC)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef WAVE_DATA
!
!-----------------------------------------------------------------------
!  Set surface wind-induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
#  ifdef ANA_WWAVE
      CALL ana_wwave (ng, tile, iADM)
#  else
#   ifdef WAVES_DIR
      CALL set_2dfldr_tile (ng, tile, iADM, idWdir,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%DwaveG,                          &
     &                      FORCES(ng)%Dwave,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

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

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Dwave)
      END IF

#     ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iADM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Dwave)
#     endif
#    endif
#   endif

#   ifdef WAVES_HEIGHT
      CALL set_2dfldr_tile (ng, tile, iADM, idWamp,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%HwaveG,                          &
     &                      FORCES(ng)%Hwave,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef WAVES_LENGTH
      CALL set_2dfldr_tile (ng, tile, iADM, idWlen,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%LwaveG,                          &
     &                      FORCES(ng)%Lwave,                           &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef WAVES_TOP_PERIOD
      CALL set_2dfldr_tile (ng, tile, iADM, idWptp,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%Pwave_topG,                      &
     &                      FORCES(ng)%Pwave_top,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef WAVES_BOT_PERIOD
      CALL set_2dfldr_tile (ng, tile, iADM, idWpbt,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%Pwave_botG,                      &
     &                      FORCES(ng)%Pwave_bot,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined WAVES_UB
      CALL set_2dfldr_tile (ng, tile, iADM, idWorb,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%Ub_swanG,                        &
     &                      FORCES(ng)%Ub_swan,                         &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined TKE_WAVEDISS
      CALL set_2dfldr_tile (ng, tile, iADM, idWdis,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%Wave_dissipG,                    &
     &                      FORCES(ng)%Wave_dissip,                     &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined SVENDSEN_ROLLER
      CALL set_2dfldr_tile (ng, tile, iADM, idWbrk,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng)%Wave_breakG,                     &
     &                      FORCES(ng)%Wave_break,                      &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
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
      CALL ana_specir (ng, tile, iADM)
# endif

# ifdef ANA_SPINNING
!
!-----------------------------------------------------------------------
!  Set time-varying rotation force (centripetal accelerations) for
!  polar coordinate grids.
!-----------------------------------------------------------------------
!
      CALL ana_spinning (ng, tile, iADM)
# endif
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
!
# ifdef ANA_PSOURCE
      IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
        CALL ana_psource (ng, tile, iADM)
      END IF
# else
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (LuvSrc(ng).or.LwSrc(ng)) THEN
          CALL set_ngfldr (ng, iADM, idRtra, 1, Nsrc(ng), 1,            &
     &                     1, Nsrc(ng), 1,                              &
     &                     SOURCES(ng) % QbarG,                         &
     &                     SOURCES(ng) % Qbar,                          &
     &                     update)
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
            CALL set_ngfldr (ng, iADM, idRtrc(itrc), 1, Nsrc(ng), N(ng),&
     &                       1, Nsrc(ng), N(ng),                        &
     &                       SOURCES(ng) % TsrcG(:,:,:,itrc),           &
     &                       SOURCES(ng) % Tsrc(:,:,itrc),              &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.  In grid refinement, only the
!  coarser grid (RefineScale(ng)=0) open boundary conditions data is
!  processed and needed.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
# ifdef ANA_FSOBC
        CALL ana_fsobc (ng, tile, iADM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idZbry(iwest), JLB, JUB, 1,      &
     &                       0, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % zetaG_west,                 &
     &                       BOUNDARY(ng) % zeta_west,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idZbry(ieast), JLB, JUB, 1,      &
     &                       0, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % zetaG_east,                 &
     &                       BOUNDARY(ng) % zeta_east,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idZbry(isouth), ILB, IUB, 1,     &
     &                       0, Lm(ng)+1 ,1,                            &
     &                       BOUNDARY(ng) % zetaG_south,                &
     &                       BOUNDARY(ng) % zeta_south,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idZbry(inorth), ILB, IUB, 1,     &
     &                       0, Lm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % zetaG_north,                &
     &                       BOUNDARY(ng) % zeta_north,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif

# if defined WET_DRY
!
!  Ensure that water level on boundary cells is above bed elevation.
!
        IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
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
          CALL mp_boundary (ng, iADM, JstrR, JendR, JLB, JUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_west)
#  endif
      END IF

        IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
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
          CALL mp_boundary (ng, iADM, JstrR, JendR, JLB, JUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_east)
#  endif
        END IF

        IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
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
          CALL mp_boundary (ng, iADM, IstrR, IendR, ILB, IUB, 1, 1,     &
     &                      bry_update,                                 &
     &                      BOUNDARY(ng)%zeta_south)
#  endif
        END IF

        IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
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
          CALL mp_boundary (ng, iADM, IstrR, IendR, ILB, IUB, 1, 1,     &
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
        CALL ana_m2obc (ng, tile, iADM)
# else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (ad_LBC(iwest,isUbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU2bc(iwest), JLB, JUB, 1,      &
     &                       0, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % ubarG_west,                 &
     &                       BOUNDARY(ng) % ubar_west,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(iwest,isVbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV2bc(iwest), JLB, JUB, 1,      &
     &                       1, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % vbarG_west,                 &
     &                       BOUNDARY(ng) % vbar_west,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(ieast,isUbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU2bc(ieast), JLB, JUB, 1,      &
     &                       0, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % ubarG_east,                 &
     &                       BOUNDARY(ng) % ubar_east,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(ieast,isVbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV2bc(ieast), JLB, JUB, 1,      &
     &                       1, Mm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % vbarG_east,                 &
     &                       BOUNDARY(ng) % vbar_east,                  &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(isouth,isUbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU2bc(isouth), ILB, IUB, 1,     &
     &                       1, Lm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % ubarG_south,                &
     &                       BOUNDARY(ng) % ubar_south,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(isouth,isVbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV2bc(isouth), ILB, IUB, 1,     &
     &                       0, Lm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % vbarG_south,                &
     &                       BOUNDARY(ng) % vbar_south,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(inorth,isUbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU2bc(inorth), ILB, IUB, 1,     &
     &                       1, Lm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % ubarG_north,                &
     &                       BOUNDARY(ng) % ubar_north,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(inorth,isVbar,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV2bc(inorth), ILB, IUB, 1,     &
     &                       0, Lm(ng)+1, 1,                            &
     &                       BOUNDARY(ng) % vbarG_north,                &
     &                       BOUNDARY(ng) % vbar_north,                 &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
# endif
      END IF

# ifdef SOLVE3D
!
!  3D momentum.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
#  ifdef ANA_M3OBC
        CALL ana_m3obc (ng, tile, iADM)
#  else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (ad_LBC(iwest,isUvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU3bc(iwest), JLB, JUB, N(ng),  &
     &                       0, Mm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % uG_west,                    &
     &                       BOUNDARY(ng) % u_west,                     &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(iwest,isVvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV3bc(iwest), JLB, JUB, N(ng),  &
     &                       1, Mm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % vG_west,                    &
     &                       BOUNDARY(ng) % v_west,                     &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(ieast,isUvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU3bc(ieast), JLB, JUB, N(ng),  &
     &                       0, Mm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % uG_east,                    &
     &                       BOUNDARY(ng) % u_east,                     &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(ieast,isVvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV3bc(ieast), JLB, JUB, N(ng),  &
     &                       1, Mm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % vG_east,                    &
     &                       BOUNDARY(ng) % v_east,                     &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(isouth,isUvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU3bc(isouth), ILB, IUB, N(ng), &
     &                       1, Lm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % uG_south,                   &
     &                       BOUNDARY(ng) % u_south,                    &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(isouth,isVvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV3bc(isouth), ILB, IUB, N(ng), &
     &                       0, Lm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % vG_south,                   &
     &                       BOUNDARY(ng) % v_south,                    &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(inorth,isUvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idU3bc(inorth), ILB, IUB, N(ng), &
     &                       1, Lm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % uG_north,                   &
     &                       BOUNDARY(ng) % u_north,                    &
     &                       update)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (ad_LBC(inorth,isUvel,ng)%acquire) THEN
            CALL set_ngfldr (ng, iADM, idV3bc(inorth), ILB, IUB, N(ng), &
     &                       0, Lm(ng)+1, N(ng),                        &
     &                       BOUNDARY(ng) % vG_north,                   &
     &                       BOUNDARY(ng) % v_north,                    &
     &                       update)
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
        CALL ana_tobc (ng, tile, iADM)
#  else
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          DO itrc=1,NT(ng)
            IF (ad_LBC(iwest,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfldr (ng, iADM, idTbry(iwest,itrc),            &
     &                         JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),     &
     &                         BOUNDARY(ng) % tG_west(:,:,:,itrc),      &
     &                         BOUNDARY(ng) % t_west(:,:,itrc),         &
     &                         update)
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (ad_LBC(ieast,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfldr (ng, iADM, idTbry(ieast,itrc),            &
     &                         JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),     &
     &                         BOUNDARY(ng) % tG_east(:,:,:,itrc),      &
     &                         BOUNDARY(ng) % t_east(:,:,itrc),         &
     &                         update)
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (ad_LBC(isouth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfldr (ng, iADM, idTbry(isouth,itrc),           &
     &                         ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),     &
     &                         BOUNDARY(ng) % tG_south(:,:,:,itrc),     &
     &                         BOUNDARY(ng) % t_south(:,:,itrc),        &
     &                         update)
              IF (exit_flag.ne.NoError) RETURN
            END IF

            IF (ad_LBC(inorth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfldr (ng, iADM, idTbry(inorth,itrc),           &
     &                         ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),     &
     &                         BOUNDARY(ng) % tG_north(:,:,:,itrc),     &
     &                         BOUNDARY(ng) % t_north(:,:,itrc),        &
     &                         update)
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
        CALL ana_ssh (ng, tile, iADM)
# else
        CALL set_2dfldr_tile (ng, tile, iADM, idSSHc,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CLIMA(ng)%sshG,                           &
     &                        CLIMA(ng)%ssh,                            &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
      IF (Lm2CLM(ng)) THEN
# ifdef ANA_M2CLIMA
        CALL ana_m2clima (ng, tile, iADM)
# else
        CALL set_2dfldr_tile (ng, tile, iADM, idUbcl,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CLIMA(ng)%ubarclmG,                       &
     &                        CLIMA(ng)%ubarclm,                        &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL set_2dfldr_tile (ng, tile, iADM, idVbcl,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CLIMA(ng)%vbarclmG,                       &
     &                        CLIMA(ng)%vbarclm,                        &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
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
        CALL ana_m3clima (ng, tile, iADM)
#  else
        CALL set_3dfldr_tile (ng, tile, iADM, idUclm,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                        CLIMA(ng)%uclmG,                            &
     &                        CLIMA(ng)%uclm,                             &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL set_3dfldr_tile (ng, tile, iADM, idVclm,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                        CLIMA(ng)%vclmG,                            &
     &                        CLIMA(ng)%vclm,                             &
     &                        update)

        IF (exit_flag.ne.NoError) RETURN
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#  ifdef ANA_TCLIMA
      IF (ANY(LtracerCLM(:,ng))) THEN
        CALL ana_tclima (ng, tile, iADM)
      END IF
#  else
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng)) THEN
          ic=ic+1
          CALL set_3dfldr_tile (ng, tile, iADM, idTclm(itrc),           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          CLIMA(ng)%tclmG(:,:,:,:,ic),            &
     &                          CLIMA(ng)%tclm (:,:,:,ic),              &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
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
        CALL set_2dfldr_tile (ng, tile, iADM, idFsur,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zetaG,                          &
     &                        OCEAN(ng)%zeta(:,:,k),                    &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
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
        CALL set_2dfldr_tile (ng, tile, iADM, idUbar,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubarG,                          &
     &                        OCEAN(ng)%ubar(:,:,k),                    &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_2dfldr_tile (ng, tile, iADM, idVbar,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbarG,                          &
     &                        OCEAN(ng)%vbar(:,:,k),                    &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO

#  ifdef FORWARD_RHS
!
!  Set forward variables associated with 2D right-hand-side terms.
!
      DO k=1,2
        CALL set_2dfldr_tile (ng, tile, iADM, idRzet,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%rzetaG,                         &
     &                        OCEAN(ng)%rzeta(:,:,k),                   &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_2dfldr_tile (ng, tile, iADM, idRu2d,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%rubarG,                         &
     &                        OCEAN(ng)%rubar(:,:,k),                   &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_2dfldr_tile (ng, tile, iADM, idRv2d,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%rvbarG,                         &
     &                        OCEAN(ng)%rvbar(:,:,k),                   &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  endif

#  ifdef SOLVE3D
!
!  Set forward time-averaged barotropic fluxes.
!
      CALL set_2dfldr_tile (ng, tile, iADM, idUfx1,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      COUPLING(ng)%DU_avg1G,                      &
     &                      COUPLING(ng)%DU_avg1,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfldr_tile (ng, tile, iADM, idUfx2,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      COUPLING(ng)%DU_avg2G,                      &
     &                      COUPLING(ng)%DU_avg2,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfldr_tile (ng, tile, iADM, idVfx1,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      COUPLING(ng)%DV_avg1G,                      &
     &                      COUPLING(ng)%DV_avg1,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

      CALL set_2dfldr_tile (ng, tile, iADM, idVfx2,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      COUPLING(ng)%DV_avg2G,                      &
     &                      COUPLING(ng)%DV_avg2,                       &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set 3D momentum.
!
      DO k=1,2
        CALL set_3dfldr_tile (ng, tile, iADM, idUvel,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%uG,                             &
     &                        OCEAN(ng)%u(:,:,:,k),                     &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_3dfldr_tile (ng, tile, iADM, idVvel,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%vG,                             &
     &                        OCEAN(ng)%v(:,:,:,k),                     &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO

#   ifdef FORWARD_RHS
!
!  Set variables associated with 3D momentum right-hand-side terms.
!
      DO k=1,2
        CALL set_3dfldr_tile (ng, tile, iADM, idRu3d,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%ruG,                            &
     &                        OCEAN(ng)%ru(:,:,:,k),                    &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN

        CALL set_3dfldr_tile (ng, tile, iADM, idRv3d,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%rvG,                            &
     &                        OCEAN(ng)%rv(:,:,:,k),                    &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif
!
!  Set 3D tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,3
          CALL set_3dfldr_tile (ng, tile, iADM, idTvar(itrc),           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          OCEAN(ng)%tG(:,:,:,:,itrc),             &
     &                          OCEAN(ng)%t(:,:,:,k,itrc),              &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END DO

#   ifdef FORWARD_MIXING
!
!  Set forward vertical mixing variables.
!
      DO itrc=1,NAT
        CALL set_3dfldr_tile (ng, tile, iADM, idDiff(itrc),             &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        MIXING(ng)%AktG(:,:,:,:,itrc),            &
     &                        MIXING(ng)%Akt(:,:,:,itrc),               &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO

      CALL set_3dfldr_tile (ng, tile, iADM, idVvis,                     &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      MIXING(ng)%AkvG,                            &
     &                      MIXING(ng)%Akv,                             &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined MY25_MIXING_NOT_YET || defined GLS_MIXING_NOT_YET
!
!  Set forward turbulent kinetic energy.
!
      DO k=1,3
        CALL set_3dfldr_tile (ng, tile, iADM, idMtke,                   &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        MIXING(ng)%tkeG,                          &
     &                        MIXING(ng)%tke(:,:,:,k),                  &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Set forward turbulent kinetic energy times length scale.
!
      DO k=1,3
        CALL set_3dfldr_tile (ng, tile, iADM, idMtls,                   &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        MIXING(ng)%glsG,                          &
     &                        MIXING(ng)%gls(:,:,:,k),                  &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Set forward vertical mixing length scale.
!
      CALL set_3dfldr_tile (ng, tile, iADM, idVmLS,                     &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      MIXING(ng)%LscaleG,                         &
     &                      MIXING(ng)%Lscale,                          &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set forward vertical mixing coefficient for turbulent kinetic energy.
!
      CALL set_3dfldr_tile (ng, tile, iADM, idVmKK,                     &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      MIXING(ng)%AkkG,                            &
     &                      MIXING(ng)%Akk,                             &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN

#    ifdef GLS_MIXING_NOT_YET
!
!  Set forward vertical mixing coefficient for turbulent length scale.
!
      CALL set_3dfldr_tile (ng, tile, iADM, idVmKP,                     &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      MIXING(ng)%AkpG,                            &
     &                      MIXING(ng)%Akp,                             &
     &                      update)

      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif

#   ifdef LMD_MIXING_NOT_YET
!
!  Set forward depth of surface oceanic boundary layer.
!
      CALL set_2dfldr_tile (ng, tile, iADM, idHsbl,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      MIXING(ng)%hsblG,                           &
     &                      MIXING(ng)%hsbl,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef LMD_BKPP_NOT_YET
!
!  Set forward depth of bottom oceanic boundary layer.
!
      CALL set_2dfldr_tile (ng, tile, iADM, idHbbl,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      MIXING(ng)%hbblG,                           &
     &                      MIXING(ng)%hbbl,                            &
     &                      update)
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef LMD_NONLOCAL_NOT_YET
!
!  Set forward boundary layer nonlocal transport.
!
      DO itrc=1,NAT
        CALL set_3dfldr_tile (ng, tile, iADM, idGhat(itrc),             &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        MIXING(ng)%ghatsG(:,:,:,:,itrc),          &
     &                        MIXING(ng)%ghat(:,:,:,itrc),              &
     &                        update)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif
#  endif
# endif

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Set forcing adjoint sensitivity functional.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
#  ifdef SENSITIVITY_4DVAR
      IF (LsenPSAS(ng)) THEN
#  endif
        IF (SCALARS(ng)%Lstate(isFsur)) THEN
          CALL set_2dfldr_tile (ng, tile, iADM, idZads,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CLIMA(ng)%zeta_adsG,                    &
     &                          CLIMA(ng)%zeta_ads,                     &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D momentum.
!
        IF (SCALARS(ng)%Lstate(isUbar)) THEN
          CALL set_2dfldr_tile (ng, tile, iADM, idUbas,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CLIMA(ng)%ubar_adsG,                    &
     &                          CLIMA(ng)%ubar_ads,                     &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
        IF (SCALARS(ng)%Lstate(isVbar)) THEN
          CALL set_2dfldr_tile (ng, tile, iADM, idVbas,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CLIMA(ng)%vbar_adsG,                    &
     &                          CLIMA(ng)%vbar_ads,                     &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Set 3D momentum.
!
        IF (SCALARS(ng)%Lstate(isUvel)) THEN
          CALL set_3dfldr_tile (ng, tile, iADM, idUads,                 &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          CLIMA(ng)%u_adsG,                       &
     &                          CLIMA(ng)%u_ads,                        &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
        IF (SCALARS(ng)%Lstate(isVvel)) THEN
          CALL set_3dfldr_tile (ng, tile, iADM, idVads,                 &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          CLIMA(ng)%v_adsG,                       &
     &                          CLIMA(ng)%v_ads,                        &
     &                          update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Set 3D tracers.
!
        DO itrc=1,NT(ng)
          IF (SCALARS(ng)%Lstate(isTvar(itrc))) THEN
            CALL set_3dfldr_tile (ng, tile, iADM, idTads(itrc),         &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng)%t_adsG(:,:,:,:,itrc),       &
     &                            CLIMA(ng)%t_ads(:,:,:,itrc),          &
     &                            update)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
#  ifdef SENSITIVITY_4DVAR
      END IF
#  endif
# endif

# if defined ADJUST_STFLUX    || defined ADJUST_WSTRESS      || \
     defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Clear adjoint vertical boundary conditions arrays for IO purposes at
!  every time-step. Otherwise, their cumulative values will be written
!  instead of instantaneous values.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          FORCES(ng)%ad_sustr(i,j)=0.0_r8
          FORCES(ng)%ad_svstr(i,j)=0.0_r8
        END DO
      END DO
#  ifdef SOLVE3D
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            FORCES(ng)%ad_stflx(i,j,itrc)=0.0_r8
            FORCES(ng)%ad_btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END DO
#   ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          FORCES(ng)%ad_srflx(i,j)=0.0_r8
        END DO
      END DO
#   endif
#   ifdef BULK_FLUXES
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          FORCES(ng)%ad_lhflx(i,j)=0.0_r8
          FORCES(ng)%ad_shflx(i,j)=0.0_r8
          FORCES(ng)%ad_lrflx(i,j)=0.0_r8
#    ifdef EMINUSP
          FORCES(ng)%ad_evap(i,j)=0.0_r8
#    endif
        END DO
      END DO
#   endif
#  endif
# endif

      RETURN
      END SUBROUTINE ad_set_data_tile
#else
      SUBROUTINE ad_set_data
      RETURN
      END SUBROUTINE ad_set_data
#endif
