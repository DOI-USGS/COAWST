#include "cppdefs.h"
      MODULE mod_forces
!
!git $Id$
!svn $Id: mod_forces.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Surface momentum stresses.                                          !
!                                                                      !
!  sustr        Surface momentum flux (wind stress) in the             !
!                 XI-direction (m2/s2) at horizontal U-points.         !
!  sustrG       Latest two-time snapshots of input "sustr" grided      !
!                 data used for interpolation.                         !
!  svstr        Surface momentum flux (wind stress) in the             !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!  svstrG       Latest two-time snapshots of input "svstr" grided      !
!                 data used for interpolation.                         !
!  Taux         Surface stress in the XI-direction at rho points       !
!                 from atm model.                                      !
!  Tauy         Surface stress in the ETA-direction at rho points      !
!                 from atm model.                                      !
!                                                                      !
!  Bottom momentum stresses.                                           !
!                                                                      !
!  bustr        Bottom momentum flux (bottom stress) in the            !
!                 XI-direction (m2/s2) at horizontal U-points.         !
!  bvstr        Bottom momentum flux (bottom stress) in the            !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!                                                                      !
!  Surface wind induced waves.                                         !
!                                                                      !
!  Hwave        Surface wind induced wave height (m).                  !
!  HwaveG       Latest two-time snapshots of input "Hwave" grided      !
!                 data used for interpolation.                         !
!  Dwave        Surface wind induced mean wave direction (radians).    !
!  DwaveG       Latest two-time snapshots of input "Dwave" grided      !
!                 data used for interpolation.                         !
!  Dwavep       Surface wind induced peak wave direction (radians).    !
!  DwavepG      Latest two-time snapshots of input "Dwavep" grided     !
!                 data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output       !
!  LwaveG       Latest two-time snapshots of input "Lwave" grided      !
!                 data used for interpolation.                         !
!  Lwavep       Peak surface wavelength read in from swan output       !
!  LwavepG      Latest two-time snapshots of input "Lwavep" grided     !
!                 data used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).                  !
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided  !
!                 data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).                   !
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided  !
!                 data used for interpolation.                         !
!  Uwave_rms    Bottom orbital velocity read in from swan output       !
!  Uwave_rmsG   Latest two-time snapshots of input "Uwave_rms" grided  !
!                 data used for interpolation.                         !
!  wave_dissip  Wave dissipation                                       !
!  wave_dissipG Latest two-time snapshots of input "wave_dissip"       !
!                 gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.    !
!  Wave_breakG  Latest two-time snapshots of input "wave_break"        !
!                 gridded data used for interpolation.                 !
!  Wave_ds      Wave directional spreading.                            !
!  Wave_qp      Wave spectrum peakedness.                              !
!                                                                      !
!  Solar shortwave radiation flux.                                     !
!                                                                      !
!  srflx        Surface shortwave solar radiation flux (degC m/s)      !
!                  at horizontal RHO-points                            !
!  srflxG       Latest two-time snapshots of input "srflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Cloud fraction.                                                     !
!                                                                      !
!  cloud        Cloud fraction (percentage/100).                       !
!  cloudG       Latest two-time snapshots of input "cloud" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat fluxes, Atmosphere-Ocean bulk parameterization.        !
!                                                                      !
!  lhflx        Latent heat flux (degC m/s).                           !
!  lrflx        Longwave radiation (degC m/s).                         !
!  shflx        Sensible heat flux (degC m/s).                         !
!                                                                      !
!  Surface air humidity.                                               !
!                                                                      !
!  Hair         Surface air specific (g/kg) or relative humidity       !
!                 (percentage).                                        !
!  HairG        Latest two-time snapshots of input "Hair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air pressure.                                               !
!                                                                      !
!  Pair         Surface air pressure (mb).                             !
!  PairG        Latest two-time snapshots of input "Pair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air temperature.                                            !
!                                                                      !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  Surface Winds.                                                      !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Rain fall rate.                                                     !
!                                                                      !
!  evap         Evaporation rate (kg/m2/s).                            !
!  rain         Rain fall rate (kg/m2/s).                              !
!  rainG        Latest two-time snapshots of input "rain" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface tracer fluxes.                                              !
!                                                                      !
!  stflux       Forcing surface flux of tracer type variables from     !
!                 data, coupling, bulk flux parameterization, or       !
!                 analytical formulas.                                 !
!                                                                      !
!                 stflux(:,:,itemp)  surface net heat flux             !
!                 stflux(:,:,isalt)  surface net freshwater flux (E-P) !
!                                                                      !
!  stfluxG      Latest two-time snapshots of input "stflux" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  stflx        ROMS state surface flux of tracer type variables       !
!                 (TracerUnits m/s) at horizontal RHO-points, as used  !
!                 in the governing equations.                          !

#if defined ADJUST_STFLUX && defined FOUR_DVAR
!                                                                      !
!  tflux        Surface tracer flux adjustments (:,:,Nfrec,2,itrc)     !
#endif
!                                                                      !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflux       Forcing bottom flux of tracer type variables from      !
!                 data or analytical formulas. Usually, the bottom     !
!                 flux of tracer is zero.                              !
!                                                                      !
!                 btflux(:,:,itemp)  bottom heat flux                  !
!                 btflux(:,:,isalt)  bottom freshwater flux            !
!                                                                      !
!  btfluxG      Latest two-time snapshots of input "vtflux" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  btflx        ROMS state bottom flux of tracer type variables        !
!                 (TracerUnits m/s) at horizontal RHO-points, as used  !
!                 in the governing equations.                          !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Surface net heat flux sensitivity to SST,              !
!                 d(Q)/d(SST), (m/s).                                  !
!  dqdtG        Latest two-time snapshots of input "dqdt" grided       !
!                 data used for interpolation.                         !
!  sst          Sea surface temperature (Celsius).                     !
!  sstG         Latest two-time snapshots of input "sst" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface freshwater flux correction.                                 !
!                                                                      !
!  sss          Sea surface salinity (PSU).                            !
!  sssG         Latest two-time snapshots of input "sss" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface spectral downwelling irradiance.                            !
!                                                                      !
!  SpecIr       Spectral irradiance (NBands) from 400-700 nm at        !
!                 5 nm bandwidth.                                      !
!  avcos        Cosine of average zenith angle of downwelling          !
!                 spectral photons.                                    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
#if defined WAVE_PARTITION
# if defined INWAVE_MODEL
        USE mod_inwave_params
# endif
# if defined SWAN_MODEL
        USE mct_coupler_params
# endif
#endif
!
        implicit none
!
        PUBLIC :: allocate_forces
        PUBLIC :: deallocate_forces
        PUBLIC :: initialize_forces
!
!-----------------------------------------------------------------------
!  Define T_FORCES structure.
!-----------------------------------------------------------------------
!
        TYPE T_FORCES
!
!  Nonlinear model state.
!
          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
#if !defined ANA_SMFLUX  && !defined BULK_FLUXES || \
     defined FORWARD_FLUXES
          real(r8), pointer :: sustrG(:,:,:)
          real(r8), pointer :: svstrG(:,:,:)
#endif
#if defined ATM2OCN_FLUXES
          real(r8), pointer :: Taux(:,:)
          real(r8), pointer :: Tauy(:,:)
#endif
#if defined WAV2OCN_FLUXES
          real(r8), pointer :: Tauocx(:,:)
          real(r8), pointer :: Tauocy(:,:)
#endif
#ifdef ADJUST_WSTRESS
          real(r8), pointer :: ustr(:,:,:,:)
          real(r8), pointer :: vstr(:,:,:,:)
#endif
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
    defined SPECTRAL_LIGHT || defined AIR_OCEAN
          real(r8), pointer :: Pair(:,:)
# ifndef ANA_PAIR
          real(r8), pointer :: PairG(:,:,:)
# endif
#endif

#ifdef WAVES_DIR
          real(r8), pointer :: Dwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: DwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_DIRP
          real(r8), pointer :: Dwavep(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: DwavepG(:,:,:)
# endif
#endif

#ifdef WAVES_HEIGHT
          real(r8), pointer :: Hwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: HwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_LENGTH
          real(r8), pointer :: Lwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: LwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_LENGTHP
          real(r8), pointer :: Lwavep(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: LwavepG(:,:,:)
# endif
#endif

#ifdef WAVES_TOP_PERIOD
          real(r8), pointer :: Pwave_top(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Pwave_topG(:,:,:)
# endif
#endif

#ifdef WAVES_BOT_PERIOD
          real(r8), pointer :: Pwave_bot(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Pwave_botG(:,:,:)
# endif
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN || \
    defined SED_BEDLOAD_VANDERA
          real(r8), pointer :: Uwave_rms(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Uwave_rmsG(:,:,:)
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
# if defined WAVE_PARTITION
          real(r8), pointer :: Dissip_break_dir(:,:,:)
# endif
          real(r8), pointer :: Dissip_break(:,:)
# if defined DISSIP_BREAK_DIR
          real(r8), pointer :: Dissip_breakx(:,:)
          real(r8), pointer :: Dissip_breaky(:,:)
          real(r8), pointer :: Dissip_wcapx(:,:)
          real(r8), pointer :: Dissip_wcapy(:,:)
# endif
          real(r8), pointer :: Dissip_wcap(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_breakG(:,:,:)
          real(r8), pointer :: Dissip_wcapG(:,:,:)
# endif
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          real(r8), pointer :: Dissip_fric(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_fricG(:,:,:)
# endif
#endif

#if defined WEC_ROLLER
          real(r8), pointer :: Dissip_roller(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_rollerG(:,:,:)
# endif
#endif

#if defined ROLLER_SVENDSEN || defined WW3_MODEL
          real(r8), pointer :: wave_break(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: wave_breakG(:,:,:)
# endif
#endif

#if defined WAVES_DSPR
          real(r8), pointer :: Wave_ds(:,:)
          real(r8), pointer :: Wave_qp(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Wave_dsG(:,:,:)
          real(r8), pointer :: Wave_qpG(:,:,:)
# endif
#endif

#if defined WEC_ROLLER
          real(r8), pointer :: rollA(:,:)
#endif

#ifdef SPECTRUM_STOKES
          real(r8), pointer :: spec_wn(:,:,:)
          real(r8), pointer :: spec_us(:,:,:)
          real(r8), pointer :: spec_vs(:,:,:)
#endif

#if defined WAVE_PARTITION
          real(r8), pointer :: HsPt01(:,:)
          real(r8), pointer :: HsPt02(:,:)
          real(r8), pointer :: TpPt01(:,:)
          real(r8), pointer :: TpPt02(:,:)
          real(r8), pointer :: WlPt01(:,:)
          real(r8), pointer :: WlPt02(:,:)
          real(r8), pointer :: DrPt01(:,:)
          real(r8), pointer :: DrPt02(:,:)
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
          real(r8), pointer :: srflx(:,:)
#  ifndef ANA_SRFLUX
          real(r8), pointer :: srflxG(:,:,:)
#  endif
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
          real(r8), pointer :: srflx_avg(:,:)
          real(r8), pointer :: srflxG_avg(:,:,:)
# endif

# ifdef CLOUDS
          real(r8), pointer :: cloud(:,:)
#  ifndef ANA_CLOUD
          real(r8), pointer :: cloudG(:,:,:)
#  endif
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lrflx(:,:)
          real(r8), pointer :: shflx(:,:)
# endif

# if !defined LONGWAVE && (defined BULK_FLUXES || \
      defined AIR_OCEAN)
          real(r8), pointer :: lrflxG(:,:,:)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined AIR_OCEAN || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
          real(r8), pointer :: Hair(:,:)
#  ifndef ANA_HUMIDITY
          real(r8), pointer :: HairG(:,:,:)
#  endif
          real(r8), pointer :: Tair(:,:)
#  ifndef ANA_TAIR
          real(r8), pointer :: TairG(:,:,:)
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
     defined AIR_OCEAN ||  defined SPECTRAL_LIGHT
          real(r8), pointer :: Uwind(:,:)
          real(r8), pointer :: Vwind(:,:)
#  ifndef ANA_WINDS
          real(r8), pointer :: UwindG(:,:,:)
          real(r8), pointer :: VwindG(:,:,:)
#  endif
# endif

# if defined BULK_FLUXES
          real(r8), pointer :: Cd(:,:)
# endif

# if defined BULK_FLUXES || defined AIR_OCEAN
          real(r8), pointer :: rain(:,:)
#  ifndef ANA_RAIN
          real(r8), pointer :: rainG(:,:,:)
#  endif
#  ifdef EMINUSP
          real(r8), pointer :: evap(:,:)
#  endif
# endif

          real(r8), pointer :: stflux(:,:,:)
# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
          real(r8), pointer :: stfluxG(:,:,:,:)
# endif
          real(r8), pointer :: stflx(:,:,:)
# ifdef ADJUST_STFLUX
          real(r8), pointer :: tflux(:,:,:,:,:)
# endif

          real(r8), pointer :: btflux(:,:,:)
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
          real(r8), pointer :: btfluxG(:,:,:,:)
# endif
          real(r8), pointer :: btflx(:,:,:)

# ifdef QCORRECTION
          real(r8), pointer :: dqdt(:,:)
          real(r8), pointer :: sst(:,:)
#  ifndef ANA_SST
          real(r8), pointer :: dqdtG(:,:,:)
          real(r8), pointer :: sstG(:,:,:)
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
          real(r8), pointer :: sss(:,:)
#  ifndef ANA_SSS
          real(r8), pointer :: sssG(:,:,:)
#  endif
# endif

# if defined ECOSIM  || defined SPECTRAL_LIGHT
          real(r8), pointer :: SpecIr(:,:,:)
          real(r8), pointer :: avcos(:,:,:)
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_sustr(:,:)
          real(r8), pointer :: tl_svstr(:,:)
# ifdef ADJUST_WSTRESS
          real(r8), pointer :: tl_ustr(:,:,:,:)
          real(r8), pointer :: tl_vstr(:,:,:,:)
# endif
          real(r8), pointer :: tl_bustr(:,:)
          real(r8), pointer :: tl_bvstr(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: tl_stflux(:,:,:)
          real(r8), pointer :: tl_stflx (:,:,:)
          real(r8), pointer :: tl_btflux(:,:,:)
          real(r8), pointer :: tl_btflx (:,:,:)
#  ifdef ADJUST_STFLUX
          real(r8), pointer :: tl_tflux(:,:,:,:,:)
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: tl_srflx(:,:)
#  endif
#  ifdef BULK_FLUXES
          real(r8), pointer :: tl_lhflx(:,:)
          real(r8), pointer :: tl_lrflx(:,:)
          real(r8), pointer :: tl_shflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: tl_evap(:,:)
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_sustr(:,:)
          real(r8), pointer :: ad_svstr(:,:)
# ifdef ADJUST_WSTRESS
          real(r8), pointer :: ad_ustr(:,:,:,:)
          real(r8), pointer :: ad_vstr(:,:,:,:)
# endif
          real(r8), pointer :: ad_bustr(:,:)
          real(r8), pointer :: ad_bvstr(:,:)
          real(r8), pointer :: ad_bustr_sol(:,:)
          real(r8), pointer :: ad_bvstr_sol(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: ad_stflx(:,:,:)
          real(r8), pointer :: ad_btflx(:,:,:)
#  ifdef ADJUST_STFLUX
          real(r8), pointer :: ad_tflux(:,:,:,:,:)
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: ad_srflx(:,:)
#  endif
#  ifdef BULK_FLUXES
          real(r8), pointer :: ad_lhflx(:,:)
          real(r8), pointer :: ad_lrflx(:,:)
          real(r8), pointer :: ad_shflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: ad_evap(:,:)
#   endif
#  endif
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
#  ifdef ADJUST_WSTRESS
          real(r8), pointer :: b_sustr(:,:)
          real(r8), pointer :: b_svstr(:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
          real(r8), pointer :: b_stflx(:,:,:)
#  endif
#  ifdef FOUR_DVAR
#   ifdef ADJUST_WSTRESS
          real(r8), pointer :: d_sustr(:,:,:)
          real(r8), pointer :: d_svstr(:,:,:)
          real(r8), pointer :: e_sustr(:,:)
          real(r8), pointer :: e_svstr(:,:)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
          real(r8), pointer :: d_stflx(:,:,:,:)
          real(r8), pointer :: e_stflx(:,:,:)
#   endif
#  endif
# endif
#endif

        END TYPE T_FORCES
!
        TYPE (T_FORCES), allocatable :: FORCES(:)
!
      CONTAINS
!
      SUBROUTINE allocate_forces (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#ifdef BIOLOGY
      USE mod_biology
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Nonlinear model state
!
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#if defined ATM2OCN_FLUXES
      allocate ( FORCES(ng) % Taux(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Tauy(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif
#if defined WAV2OCN_FLUXES
      allocate ( FORCES(ng) % Tauocx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Tauocy(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif
#if !defined ANA_SMFLUX     && !defined BULK_FLUXES || \
     defined FORWARD_FLUXES
      allocate ( FORCES(ng) % sustrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( FORCES(ng) % svstrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#endif
#ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
#endif
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
    defined SPECTRAL_LIGHT || defined AIR_OCEAN
      allocate ( FORCES(ng) % Pair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_PAIR
      allocate ( FORCES(ng) % PairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_DIR
      allocate ( FORCES(ng) % Dwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % DwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_DIRP
      allocate ( FORCES(ng) % Dwavep(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % DwavepG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_HEIGHT
      allocate ( FORCES(ng) % Hwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % HwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_LENGTH
      allocate ( FORCES(ng) % Lwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % LwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_LENGTHP
      allocate ( FORCES(ng) % Lwavep(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % LwavepG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_TOP_PERIOD
      allocate ( FORCES(ng) % Pwave_top(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_topG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_BOT_PERIOD
      allocate ( FORCES(ng) % Pwave_bot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_botG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN
      allocate ( FORCES(ng) % Uwave_rms(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Uwave_rmsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
# if defined WAVE_PARTITION
      allocate ( FORCES(ng) % Dissip_break_dir(LBi:UBi,LBj:UBj,ND) )
      Dmem(ng)=Dmem(ng)+size2d*ND
# endif
      allocate ( FORCES(ng) % Dissip_break(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# if defined DISSIP_BREAK_DIR
      allocate ( FORCES(ng) % Dissip_breakx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Dissip_breaky(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Dissip_wcapx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Dissip_wcapy(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif
      allocate ( FORCES(ng) % Dissip_wcap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_breakG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % Dissip_wcapG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      allocate ( FORCES(ng) % Dissip_fric(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_fricG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WEC_ROLLER
      allocate ( FORCES(ng) % Dissip_roller(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_rollerG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined ROLLER_SVENDSEN || defined WW3_MODEL
      allocate ( FORCES(ng) % wave_break(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_breakG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WAVES_DSPR
      allocate ( FORCES(ng) % Wave_ds(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Wave_qp(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_dsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % Wave_qpG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WEC_ROLLER
      allocate ( FORCES(ng) % rollA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif

#ifdef SPECTRUM_STOKES
      allocate ( FORCES(ng) % spec_wn(LBi:UBi,LBj:UBj,MSCs) )
      Dmem(ng)=Dmem(ng)+size2d*REAL(MSCs,r8)
      allocate ( FORCES(ng) % spec_us(LBi:UBi,LBj:UBj,MSCs) )
      Dmem(ng)=Dmem(ng)+size2d*REAL(MSCs,r8)
      allocate ( FORCES(ng) % spec_vs(LBi:UBi,LBj:UBj,MSCs) )
      Dmem(ng)=Dmem(ng)+size2d*REAL(MSCs,r8)
#endif

#if defined WAVE_PARTITION
      allocate ( FORCES(ng) % HsPt01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % HsPt02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % TpPt01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % TpPt02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % WlPt01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % WlPt02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % DrPt01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % DrPt02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SRFLUX
      allocate ( FORCES(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
      allocate ( FORCES(ng) % srflx_avg(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % srflxG_avg(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif

# ifdef CLOUDS
      allocate ( FORCES(ng) % cloud(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_CLOUD
      allocate ( FORCES(ng) % cloudG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined BULK_FLUXES || defined AIR_OCEAN
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# if !defined LONGWAVE && (defined BULK_FLUXES || \
      defined AIR_OCEAN)
      allocate ( FORCES(ng) % lrflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif

# if defined BULK_FLUXES || defined ECOSIM || defined AIR_OCEAN || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
      allocate ( FORCES(ng) % Hair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_HUMIDITY
      allocate ( FORCES(ng) % HairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif

      allocate ( FORCES(ng) % Tair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_TAIR
      allocate ( FORCES(ng) % TairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
     defined AIR_OCEAN ||  defined SPECTRAL_LIGHT
      allocate ( FORCES(ng) % Uwind(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % Vwind(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_WINDS
      allocate ( FORCES(ng) % UwindG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( FORCES(ng) % VwindG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined BULK_FLUXES
      allocate ( FORCES(ng) % Cd(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# if defined BULK_FLUXES || defined AIR_OCEAN
      allocate ( FORCES(ng) % rain(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_RAIN
      allocate ( FORCES(ng) % rainG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

#  endif
#  ifdef EMINUSP
      allocate ( FORCES(ng) % evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
# endif

      allocate ( FORCES(ng) % stflux(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
      allocate ( FORCES(ng) % stfluxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NT(ng),r8)*size2d

# endif

      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

# ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tflux(LBi:UBi,LBj:UBj,nfrec(ng),          &
     &                              2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(nfrec(ng)*NT(ng),r8)*size2d

# endif

      allocate ( FORCES(ng) % btflux(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
      allocate ( FORCES(ng) % btfluxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NT(ng),r8)*size2d
# endif

      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

# ifdef QCORRECTION
      allocate ( FORCES(ng) % dqdt(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % sst(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SST
      allocate ( FORCES(ng) % dqdtG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( FORCES(ng) % sstG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      allocate ( FORCES(ng) % sss(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SSS
      allocate ( FORCES(ng) % sssG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif


# if defined ECOSIM || defined SPECTRAL_LIGHT 
      allocate ( FORCES(ng) % SpecIr(LBi:UBi,LBj:UBj,NBands) )
      Dmem(ng)=Dmem(ng)+REAL(NBands,r8)*size2d

      allocate ( FORCES(ng) % avcos(LBi:UBi,LBj:UBj,NBands) )
      Dmem(ng)=Dmem(ng)+REAL(NBands,r8)*size2d
# endif

#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state
!
      allocate ( FORCES(ng) % tl_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
# endif
      allocate ( FORCES(ng) % tl_bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef SOLVE3D
      allocate ( FORCES(ng) % tl_stflux(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_btflux(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tl_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng)*NT(ng),r8)*size2d
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % tl_srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % tl_lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   ifdef EMINUSP
      allocate ( FORCES(ng) % tl_evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state
!
      allocate ( FORCES(ng) % ad_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
# endif
      allocate ( FORCES(ng) % ad_bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bustr_sol(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bvstr_sol(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef SOLVE3D
      allocate ( FORCES(ng) % ad_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % ad_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % ad_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng)*NT(ng),r8)*size2d
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % ad_srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % ad_lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   ifdef EMINUSP
      allocate ( FORCES(ng) % ad_evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
#  endif
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
#  ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % b_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % b_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % b_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  endif
# endif
# ifdef FOUR_DVAR
#  ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % e_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % e_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % d_stflx(LBi:UBi,LBj:UBj,                  &
     &                                Nfrec(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng)*NT(ng),r8)*size2d

      allocate ( FORCES(ng) % e_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  endif
# endif
#endif
!
      RETURN
      END SUBROUTINE allocate_forces
!
      SUBROUTINE deallocate_forces (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Ngrids
#ifdef SUBOBJECT_DEALLOCATION
      USE destroy_mod, ONLY : destroy
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_forces"

#ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_FORCES structure
!  separately.
!-----------------------------------------------------------------------
!
!  Nonlinear model state
!
      IF (.not.destroy(ng, FORCES(ng)%sustr, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%svstr, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%svstr')) RETURN

# if !defined ANA_SMFLUX     && !defined BULK_FLUXES || \
      defined FORWARD_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%sustrG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%sustrG')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%svstrG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%svstrG')) RETURN
# endif

# if defined ATM2OCN_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%Taux, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%Taux')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%Tauy, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%Tauy')) RETURN
# endif

#if defined WAV2OCN_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%Tauocx, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%Taux')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Tauocy, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%Taux')) RETURN
#endif

# ifdef ADJUST_WSTRESS
      IF (.not.destroy(ng, FORCES(ng)%ustr, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%ustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%vstr, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%vstr')) RETURN
# endif

      IF (.not.destroy(ng, FORCES(ng)%bustr, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%bustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%bvstr, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%bvstr')) RETURN

#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
     defined SPECTRAL_LIGHT  || defined AIR_OCEAN
      IF (.not.destroy(ng, FORCES(ng)%Pair, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%Pair')) RETURN
#  ifndef ANA_PAIR
      IF (.not.destroy(ng, FORCES(ng)%PairG, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%PairG')) RETURN
#  endif
# endif

# ifdef WAVES_DIR
      IF (.not.destroy(ng, FORCES(ng)%Dwave, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%Dwave')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%DwaveG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%DwaveG')) RETURN
#  endif
# endif

# ifdef WAVES_DIRP
      IF (.not.destroy(ng, FORCES(ng)%Dwavep, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%Dwavep')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%DwavepG, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%DwavepG')) RETURN
#  endif
# endif

# ifdef WAVES_HEIGHT
      IF (.not.destroy(ng, FORCES(ng)%Hwave, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%Hwave')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%HwaveG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%HwaveG')) RETURN
#  endif
# endif

# ifdef WAVES_LENGTH
      IF (.not.destroy(ng, FORCES(ng)%Lwave, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%Lwave')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%LwaveG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%LwaveG')) RETURN
#  endif
# endif

# ifdef WAVES_LENGTHP
      IF (.not.destroy(ng, FORCES(ng)%Lwavep, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%Lwavep')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%LwavepG, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%LwavepG')) RETURN
#  endif
# endif

# ifdef WAVES_TOP_PERIOD
      IF (.not.destroy(ng, FORCES(ng)%Pwave_top, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%Pwave_top')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Pwave_topG, MyFile,               &
     &                 __LINE__, 'FORCES(ng)%Pwave_topG')) RETURN
#  endif
# endif

# ifdef WAVES_BOT_PERIOD
      IF (.not.destroy(ng, FORCES(ng)%Pwave_bot, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%Pwave_bot')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Pwave_botG, MyFile,               &
     &                 __LINE__, 'FORCES(ng)%Pwave_botG')) RETURN
#  endif
# endif

# if defined BBL_MODEL || defined WAVES_OCEAN || \
    defined SED_BEDLOAD_VANDERA
      IF (.not.destroy(ng, FORCES(ng)%Uwave_rms, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%Uwave_rms')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Uwave_rmsG, MyFile,               &
     &                 __LINE__, 'FORCES(ng)%Uwave_rmsG')) RETURN
#  endif
# endif

# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
      IF (.not.destroy(ng, FORCES(ng)%Dissip_break, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%Dissip_break')) RETURN
#  if defined DISSIP_BREAK_DIR
      IF (.not.destroy(ng, FORCES(ng)%Dissip_breakx, MyFile,            &
     &                 __LINE__, 'FORCES(ng)%Dissip_breakx')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Dissip_breaky, MyFile,            &
     &                 __LINE__, 'FORCES(ng)%Dissip_breaky')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Dissip_wcapx, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%Dissip_wcapx')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Dissip_wcapy, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%Dissip_wcapy')) RETURN
#  endif
      IF (.not.destroy(ng, FORCES(ng)%Dissip_wcap, MyFile,              &
     &                 __LINE__, 'FORCES(ng)%Dissip_wcap')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Dissip_breakG, MyFile,            &
     &                 __LINE__, 'FORCES(ng)%Dissip_breakG')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Dissip_wcapG, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%Dissip_wcapG')) RETURN
#  endif
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      IF (.not.destroy(ng, FORCES(ng)%Dissip_fric, MyFile,              &
     &                 __LINE__, 'FORCES(ng)%Dissip_fric')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Dissip_fricG, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%Dissip_fricG')) RETURN
#  endif
# endif

# if defined WEC_ROLLER
      IF (.not.destroy(ng, FORCES(ng)%Dissip_roller, MyFile,            &
     &                 __LINE__, 'FORCES(ng)%Dissip_roller')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%rollA, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%rollA')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Dissip_rollerG, MyFile,           &
     &                 __LINE__, 'FORCES(ng)%Dissip_rollerG')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%rollAG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%rollAG')) RETURN
#  endif
# endif

#if defined ROLLER_SVENDSEN || defined WW3_MODEL
      IF (.not.destroy(ng, FORCES(ng)%wave_break, MyFile,               &
     &                 __LINE__, 'FORCES(ng)%wave_break')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Wave_breakG, MyFile,              &
     &                 __LINE__, 'FORCES(ng)%Wave_breakG')) RETURN
#  endif
# endif

# if defined WAVES_DSPR
      IF (.not.destroy(ng, FORCES(ng)%Wave_ds, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%Wave_ds')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Wave_qp, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%Wave_qp')) RETURN
#  ifndef ANA_WWAVE
      IF (.not.destroy(ng, FORCES(ng)%Wave_dsG, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%Wwave_dsG')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%Wave_qpG, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%Wave_qpG')) RETURN
#  endif
# endif

# ifdef SPECTRUM_STOKES
      IF (.not.destroy(ng, FORCES(ng)%spec_wn, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%spec_wn')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%spec_us, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%spec_us')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%spec_vs, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%spec_vs')) RETURN
# endif

# ifdef WAVE_PARTITION
      IF (.not.destroy(ng, FORCES(ng)%HsPt01, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%HsPt01')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%HsPt02, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%HsPt02')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%TpPt01, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%TpPt01')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%TpPt02, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%TpPt02')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%WlPt01, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%WlPt01')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%WlPt02, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%WlPt02')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%DrPt01, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%DrPt01')) RETURN
      IF (.not.destroy(ng, FORCES(ng)%DrPt02, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%DrPt02')) RETURN
# endif

# ifdef SOLVE3D

#  ifdef SHORTWAVE
      IF (.not.destroy(ng, FORCES(ng)%srflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%srflx')) RETURN
#   ifndef ANA_SRFLUX
      IF (.not.destroy(ng, FORCES(ng)%srflxG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%srflxG')) RETURN
#   endif
#  endif

#  if defined RED_TIDE && defined DAILY_SHORTWAVE
      IF (.not.destroy(ng, FORCES(ng)%srflx_avg, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%srflx_avg')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%srflxG_avg, MyFile,               &
     &                 __LINE__, 'FORCES(ng)%srflxG_avg')) RETURN
#  endif

#  ifdef CLOUDS
      IF (.not.destroy(ng, FORCES(ng)%cloud, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%cloud')) RETURN
#   ifndef ANA_CLOUD
      IF (.not.destroy(ng, FORCES(ng)%cloudG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%cloudG')) RETURN
#   endif
#  endif

#  if defined BULK_FLUXES || defined AIR_OCEAN
      IF (.not.destroy(ng, FORCES(ng)%lhflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%lhflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%lrflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%lrflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%shflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%shflx')) RETURN
#  endif

# if !defined LONGWAVE && (defined BULK_FLUXES || \
      defined AIR_OCEAN)
      IF (.not.destroy(ng, FORCES(ng)%lrflxG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%lrflxG')) RETURN
#  endif

# if defined BULK_FLUXES || defined ECOSIM || defined AIR_OCEAN || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
      IF (.not.destroy(ng, FORCES(ng)%Hair, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%Hair')) RETURN

#   ifndef ANA_HUMIDITY
      IF (.not.destroy(ng, FORCES(ng)%HairG, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%HairG')) RETURN
#   endif

      IF (.not.destroy(ng, FORCES(ng)%Tair, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%Tair')) RETURN

#   ifndef ANA_TAIR
      IF (.not.destroy(ng, FORCES(ng)%TairG, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%TairG')) RETURN
#   endif
#  endif

# if defined BULK_FLUXES || defined ECOSIM || \
     defined AIR_OCEAN ||  defined SPECTRAL_LIGHT
      IF (.not.destroy(ng, FORCES(ng)%Uwind, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%Uwind')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%Vwind, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%Vwind')) RETURN

#   ifndef ANA_WINDS
      IF (.not.destroy(ng, FORCES(ng)%UwindG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%UwindG')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%VwindG, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%VwindG')) RETURN
#   endif
#  endif

# if defined BULK_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%Cd, MyFile,                       &
     &                 __LINE__, 'FORCES(ng)%Cd')) RETURN
#  endif

# if defined BULK_FLUXES || defined AIR_OCEAN
      IF (.not.destroy(ng, FORCES(ng)%rain, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%rain')) RETURN

#   ifndef ANA_RAIN
      IF (.not.destroy(ng, FORCES(ng)%rainG, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%rainG')) RETURN
#   endif

#   ifdef EMINUSP
      IF (.not.destroy(ng, FORCES(ng)%evap, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%evap')) RETURN
#   endif
#  endif

      IF (.not.destroy(ng, FORCES(ng)%stflux, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%stflux')) RETURN

#  if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
      !defined ANA_SPFLUX
      IF (.not.destroy(ng, FORCES(ng)%stfluxG, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%stfluxG')) RETURN
#  endif

      IF (.not.destroy(ng, FORCES(ng)%stflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%stflx')) RETURN

#  ifdef ADJUST_STFLUX
      IF (.not.destroy(ng, FORCES(ng)%tflux, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%tflux')) RETURN
#  endif

      IF (.not.destroy(ng, FORCES(ng)%btflux, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%btflux')) RETURN

#  if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
      !defined ANA_BPFLUX
      IF (.not.destroy(ng, FORCES(ng)%btfluxG, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%btfluxG')) RETURN
#  endif

      IF (.not.destroy(ng, FORCES(ng)%btflx, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%btflx')) RETURN

#  ifdef QCORRECTION
      IF (.not.destroy(ng, FORCES(ng)%dqdt, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%dqdt')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%sst, MyFile,                      &
     &                 __LINE__, 'FORCES(ng)%sst')) RETURN

#   ifndef ANA_SST
      IF (.not.destroy(ng, FORCES(ng)%dqdtG, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%dqdtG')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%sstG, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%sstG')) RETURN
#   endif
#  endif

#  if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      IF (.not.destroy(ng, FORCES(ng)%sss, MyFile,                      &
     &                 __LINE__, 'FORCES(ng)%sss')) RETURN

#   ifndef ANA_SSS
      IF (.not.destroy(ng, FORCES(ng)%sssG, MyFile,                     &
     &                 __LINE__, 'FORCES(ng)%sssG')) RETURN
#   endif
#  endif

# if defined ECOSIM  || defined SPECTRAL_LIGHT
      IF (.not.destroy(ng, FORCES(ng)%SpecIr, MyFile,                   &
     &                 __LINE__, 'FORCES(ng)%SpecIr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%avcos, MyFile,                    &
     &                 __LINE__, 'FORCES(ng)%avcos')) RETURN
#  endif

# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state
!
      IF (.not.destroy(ng, FORCES(ng)%tl_sustr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_svstr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_svstr')) RETURN

#  ifdef ADJUST_WSTRESS
      IF (.not.destroy(ng, FORCES(ng)%tl_ustr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%tl_ustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_vstr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%tl_vstr')) RETURN
#  endif

      IF (.not.destroy(ng, FORCES(ng)%tl_bustr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_bustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_bvstr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_bvstr')) RETURN

#  ifdef SOLVE3D
      IF (.not.destroy(ng, FORCES(ng)%tl_stflux, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%tl_stflux')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_stflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_stflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_btflux, MyFile,                &
     &                 __LINE__, 'FORCES(ng)%tl_btflux')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_btflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_btflx')) RETURN

#   ifdef ADJUST_STFLUX
      IF (.not.destroy(ng, FORCES(ng)%tl_tflux, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_tflux')) RETURN
#   endif

#   ifdef SHORTWAVE
      IF (.not.destroy(ng, FORCES(ng)%tl_srflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_srflx')) RETURN
#   endif

#   ifdef BULK_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%tl_lhflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_lhflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_lrflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_lrflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%tl_shflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%tl_shflx')) RETURN

#    ifdef EMINUSP
      IF (.not.destroy(ng, FORCES(ng)%tl_evap, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%tl_evap')) RETURN
#    endif
#   endif
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state
!
      IF (.not.destroy(ng, FORCES(ng)%ad_sustr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_svstr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_svstr')) RETURN

#  ifdef ADJUST_WSTRESS
      IF (.not.destroy(ng, FORCES(ng)%ad_ustr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%ad_ustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_vstr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%ad_vstr')) RETURN
#  endif

      IF (.not.destroy(ng, FORCES(ng)%ad_bustr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_bustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_bvstr, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_bvstr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_bustr_sol, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%ad_bustr_sol')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_bvstr_sol, MyFile,             &
     &                 __LINE__, 'FORCES(ng)%ad_bvstr_sol')) RETURN

#  ifdef SOLVE3D
      IF (.not.destroy(ng, FORCES(ng)%ad_stflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_stflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_btflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_btflx')) RETURN

#   ifdef ADJUST_STFLUX
      IF (.not.destroy(ng, FORCES(ng)%ad_tflux, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_tflux')) RETURN
#   endif

#   ifdef SHORTWAVE
      IF (.not.destroy(ng, FORCES(ng)%ad_srflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_srflx')) RETURN
#   endif

#   ifdef BULK_FLUXES
      IF (.not.destroy(ng, FORCES(ng)%ad_lhflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_lhflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_lrflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_lrflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%ad_shflx, MyFile,                 &
     &                 __LINE__, 'FORCES(ng)%ad_shflx')) RETURN

#    ifdef EMINUSP
      IF (.not.destroy(ng, FORCES(ng)%ad_evap, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%ad_evap')) RETURN
#    endif
#   endif
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
#  if defined FOUR_DVAR || defined IMPULSE
#   ifdef ADJUST_WSTRESS
      IF (.not.destroy(ng, FORCES(ng)%b_sustr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%b_sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%b_svstr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%b_svstr')) RETURN
#   endif

#   if defined ADJUST_STFLUX && defined SOLVE3D
      IF (.not.destroy(ng, FORCES(ng)%b_stflx, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%b_stflx')) RETURN
#   endif
#  endif

#  ifdef FOUR_DVAR
#   ifdef ADJUST_WSTRESS
      IF (.not.destroy(ng, FORCES(ng)%d_sustr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%d_sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%d_svstr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%d_svstr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%e_sustr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%e_sustr')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%e_svstr, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%e_svstr')) RETURN
#   endif

#   if defined ADJUST_STFLUX && defined SOLVE3D
      IF (.not.destroy(ng, FORCES(ng)%d_stflx, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%d_stflx')) RETURN

      IF (.not.destroy(ng, FORCES(ng)%e_stflx, MyFile,                  &
     &                 __LINE__, 'FORCES(ng)%e_stflx')) RETURN
#   endif
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type FORCES structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(FORCES)) deallocate ( FORCES )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_forces
!
      SUBROUTINE initialize_forces (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
#ifdef BIOLOGY
      USE mod_biology
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
#ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % ustr(i,j,k,1) = IniVal
              FORCES(ng) % ustr(i,j,k,2) = IniVal
              FORCES(ng) % vstr(i,j,k,1) = IniVal
              FORCES(ng) % vstr(i,j,k,2) = IniVal
            END DO
#endif
            FORCES(ng) % sustr(i,j) = IniVal
            FORCES(ng) % svstr(i,j) = IniVal
#if defined ATM2OCN_FLUXES
            FORCES(ng) % Taux(i,j) = IniVal
            FORCES(ng) % Tauy(i,j) = IniVal
#endif
#if defined WAV2OCN_FLUXES
            FORCES(ng) % Tauocx(i,j) = IniVal
            FORCES(ng) % Tauocy(i,j) = IniVal
#endif
#if !defined ANA_SMFLUX     && !defined BULK_FLUXES || \
     defined FORWARD_FLUXES
            FORCES(ng) % sustrG(i,j,1) = IniVal
            FORCES(ng) % sustrG(i,j,2) = IniVal
            FORCES(ng) % svstrG(i,j,1) = IniVal
            FORCES(ng) % svstrG(i,j,2) = IniVal
#endif
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal

#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
    defined SPECTRAL_LIGHT || defined AIR_OCEAN
            FORCES(ng) % Pair(i,j) = IniVal
# ifndef ANA_PAIR
            FORCES(ng) % PairG(i,j,1) = IniVal
            FORCES(ng) % PairG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_DIR
            FORCES(ng) % Dwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % DwaveG(i,j,1) = IniVal
            FORCES(ng) % DwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_DIRP
            FORCES(ng) % Dwavep(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % DwavepG(i,j,1) = IniVal
            FORCES(ng) % DwavepG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_HEIGHT
            FORCES(ng) % Hwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % HwaveG(i,j,1) = IniVal
            FORCES(ng) % HwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_LENGTH
            FORCES(ng) % Lwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % LwaveG(i,j,1) = IniVal
            FORCES(ng) % LwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_LENGTHP
            FORCES(ng) % Lwavep(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % LwavepG(i,j,1) = IniVal
            FORCES(ng) % LwavepG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_TOP_PERIOD
            FORCES(ng) % Pwave_top(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Pwave_topG(i,j,1) = IniVal
            FORCES(ng) % Pwave_topG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_BOT_PERIOD
            FORCES(ng) % Pwave_bot(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Pwave_botG(i,j,1) = IniVal
            FORCES(ng) % Pwave_botG(i,j,2) = IniVal
# endif
#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
            FORCES(ng) % Uwave_rms(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Uwave_rmsG(i,j,1) = IniVal
            FORCES(ng) % Uwave_rmsG(i,j,2) = IniVal
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
# if defined WAVE_PARTITION
            DO k=1,ND
              FORCES(ng) % Dissip_break_dir(i,j,k) = IniVal
            END DO
# endif
            FORCES(ng) % Dissip_break(i,j) = IniVal
# if defined DISSIP_BREAK_DIR
            FORCES(ng) % Dissip_breakx(i,j) = IniVal
            FORCES(ng) % Dissip_breaky(i,j) = IniVal
            FORCES(ng) % Dissip_wcapx(i,j) = IniVal
            FORCES(ng) % Dissip_wcapy(i,j) = IniVal
# endif
            FORCES(ng) % Dissip_wcap(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_breakG(i,j,1) = IniVal
            FORCES(ng) % Dissip_breakG(i,j,2) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,1) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,2) = IniVal
# endif
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
            FORCES(ng) % Dissip_fric(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_fricG(i,j,1) = IniVal
            FORCES(ng) % Dissip_fricG(i,j,2) = IniVal
# endif
#endif

#if defined WEC_ROLLER
            FORCES(ng) % Dissip_roller(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_rollerG(i,j,1) = IniVal
            FORCES(ng) % Dissip_rollerG(i,j,2) = IniVal
# endif
#endif

#if defined ROLLER_SVENDSEN || defined WW3_MODEL
            FORCES(ng) % Wave_break(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_breakG(i,j,1) = IniVal
            FORCES(ng) % Wave_breakG(i,j,2) = IniVal
# endif
#endif

#if defined WAVES_DSPR
            FORCES(ng) % Wave_ds(i,j) = IniVal
            FORCES(ng) % Wave_qp(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_dsG(i,j,1) = IniVal
            FORCES(ng) % Wave_dsG(i,j,2) = IniVal
            FORCES(ng) % Wave_qpG(i,j,1) = IniVal
            FORCES(ng) % Wave_qpG(i,j,2) = IniVal
# endif
#endif

#if defined WEC_ROLLER
            FORCES(ng) % rollA(i,j) = IniVal
#endif

#ifdef SPECTRUM_STOKES
            DO k=1,MSCs
              FORCES(ng) % spec_wn(i,j,k) = IniVal
              FORCES(ng) % spec_us(i,j,k) = IniVal
              FORCES(ng) % spec_vs(i,j,k) = IniVal
            END DO
#endif

#ifdef WAVE_PARTITION
            FORCES(ng) % HsPt01(i,j) = IniVal
            FORCES(ng) % HsPt02(i,j) = IniVal
            FORCES(ng) % TpPt01(i,j) = IniVal
            FORCES(ng) % TpPt02(i,j) = IniVal
            FORCES(ng) % WlPt01(i,j) = IniVal
            FORCES(ng) % WlPt02(i,j) = IniVal
            FORCES(ng) % DrPt01(i,j) = IniVal
            FORCES(ng) % DrPt02(i,j) = IniVal
#endif

#ifdef SOLVE3D
# ifdef SHORTWAVE
            FORCES(ng) % srflx(i,j) = IniVal
#  ifndef ANA_SRFLUX
            FORCES(ng) % srflxG(i,j,1) = IniVal
            FORCES(ng) % srflxG(i,j,2) = IniVal
#  endif
# endif
# if defined RED_TIDE && defined DAILY_SHORTWAVE
            FORCES(ng) % srflx_avg(i,j) = IniVal
            FORCES(ng) % srflxG_avg(i,j,1) = IniVal
            FORCES(ng) % srflxG_avg(i,j,2) = IniVal
# endif
# ifdef CLOUDS
            FORCES(ng) % cloud(i,j) = IniVal
#  ifndef ANA_CLOUD
            FORCES(ng) % cloudG(i,j,1) = IniVal
            FORCES(ng) % cloudG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN
            FORCES(ng) % lhflx(i,j) = IniVal
            FORCES(ng) % lrflx(i,j) = IniVal
            FORCES(ng) % shflx(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined AIR_OCEAN || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
            FORCES(ng) % Hair(i,j) = IniVal
            FORCES(ng) % Tair(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    defined SPECTRAL_LIGHT || defined AIR_OCEAN
            FORCES(ng) % Uwind(i,j) = IniVal
            FORCES(ng) % Vwind(i,j) = IniVal
# endif
# if defined BULK_FLUXES
            FORCES(ng) % Cd(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN
            FORCES(ng) % rain(i,j) = IniVal
#  ifdef EMINUSP
            FORCES(ng) % evap(i,j) = IniVal
#  endif
# endif
# if !defined LONGWAVE && (defined BULK_FLUXES || \
      defined AIR_OCEAN)
            FORCES(ng) % lrflxG(i,j,1) = IniVal
            FORCES(ng) % lrflxG(i,j,2) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
#  ifndef ANA_HUMIDITY
            FORCES(ng) % HairG(i,j,1) = IniVal
            FORCES(ng) % HairG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined AIR_OCEAN || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
#  ifndef ANA_TAIR
            FORCES(ng) % TairG(i,j,1) = IniVal
            FORCES(ng) % TairG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    defined SPECTRAL_LIGHT || defined AIR_OCEAN
#  ifndef ANA_WINDS
            FORCES(ng) % UwindG(i,j,1) = IniVal
            FORCES(ng) % UwindG(i,j,2) = IniVal
            FORCES(ng) % VwindG(i,j,1) = IniVal
            FORCES(ng) % VwindG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN
#  ifndef ANA_RAIN
            FORCES(ng) % rainG(i,j,1) = IniVal
            FORCES(ng) % rainG(i,j,2) = IniVal
#  endif
# endif
# ifdef QCORRECTION
            FORCES(ng) % dqdt(i,j) = IniVal
            FORCES(ng) % sst(i,j) = IniVal
#  ifndef ANA_SST
            FORCES(ng) % dqdtG(i,j,1) = IniVal
            FORCES(ng) % dqdtG(i,j,2) = IniVal
            FORCES(ng) % sstG(i,j,1) = IniVal
            FORCES(ng) % sstG(i,j,2) = IniVal
#  endif
# endif
# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
            FORCES(ng) % sss(i,j) = IniVal
#  ifndef ANA_SSS
            FORCES(ng) % sssG(i,j,1) = IniVal
            FORCES(ng) % sssG(i,j,2) = IniVal
#  endif
# endif
            DO itrc=1,NT(ng)
              FORCES(ng) % stflux(i,j,itrc) = IniVal
# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
              FORCES(ng) % stfluxG(i,j,1,itrc) = IniVal
              FORCES(ng) % stfluxG(i,j,2,itrc) = IniVal
# endif
              FORCES(ng) % stflx(i,j,itrc) = IniVal
# ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % tflux(i,j,k,2,itrc) = IniVal
              END DO
# endif

              FORCES(ng) % btflux(i,j,itrc) = IniVal
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
              FORCES(ng) % btfluxG(i,j,1,itrc) = IniVal
              FORCES(ng) % btfluxG(i,j,2,itrc) = IniVal
# endif
              FORCES(ng) % btflx(i,j,itrc) = IniVal
            END DO
# if defined ECOSIM || defined SPECTRAL_LIGHT 
            DO itrc=1,NBands
              FORCES(ng) % SpecIr(i,j,itrc) = IniVal
              FORCES(ng) % avcos(i,j,itrc) = IniVal
            END DO
# endif
#endif
          END DO
        END DO
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % tl_ustr(i,j,k,1) = IniVal
              FORCES(ng) % tl_ustr(i,j,k,2) = IniVal
              FORCES(ng) % tl_vstr(i,j,k,1) = IniVal
              FORCES(ng) % tl_vstr(i,j,k,2) = IniVal
            END DO
# endif
            FORCES(ng) % tl_sustr(i,j) = IniVal
            FORCES(ng) % tl_svstr(i,j) = IniVal
            FORCES(ng) % tl_bustr(i,j) = IniVal
            FORCES(ng) % tl_bvstr(i,j) = IniVal
          END DO
# ifdef SOLVE3D
#  ifdef SHORTWAVE
          DO i=Imin,Imax
            FORCES(ng) % tl_srflx(i,j) = IniVal
          END DO
#  endif
#  ifdef BULK_FLUXES
          DO i=Imin,Imax
            FORCES(ng) % tl_lhflx(i,j) = IniVal
            FORCES(ng) % tl_lrflx(i,j) = IniVal
            FORCES(ng) % tl_shflx(i,j) = IniVal
#   ifdef EMINUSP
            FORCES(ng) % tl_evap(i,j) = IniVal
#   endif
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
#  ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % tl_tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % tl_tflux(i,j,k,2,itrc) = IniVal
              END DO
#  endif
              FORCES(ng) % tl_stflux(i,j,itrc) = IniVal
              FORCES(ng) % tl_stflx (i,j,itrc) = IniVal
              FORCES(ng) % tl_btflux(i,j,itrc) = IniVal
              FORCES(ng) % tl_btflx (i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % ad_ustr(i,j,k,1) = IniVal
              FORCES(ng) % ad_ustr(i,j,k,2) = IniVal
              FORCES(ng) % ad_vstr(i,j,k,1) = IniVal
              FORCES(ng) % ad_vstr(i,j,k,2) = IniVal
            END DO
# endif
            FORCES(ng) % ad_sustr(i,j) = IniVal
            FORCES(ng) % ad_svstr(i,j) = IniVal
            FORCES(ng) % ad_bustr(i,j) = IniVal
            FORCES(ng) % ad_bvstr(i,j) = IniVal
            FORCES(ng) % ad_bustr_sol(i,j) = IniVal
            FORCES(ng) % ad_bvstr_sol(i,j) = IniVal
          END DO
# ifdef SOLVE3D
#  ifdef SHORTWAVE
          DO i=Imin,Imax
            FORCES(ng) % ad_srflx(i,j) = IniVal
          END DO
#  endif
#  ifdef BULK_FLUXES
          DO i=Imin,Imax
            FORCES(ng) % ad_lhflx(i,j) = IniVal
            FORCES(ng) % ad_lrflx(i,j) = IniVal
            FORCES(ng) % ad_shflx(i,j) = IniVal
#   ifdef EMINUSP
            FORCES(ng) % ad_evap(i,j) = IniVal
#   endif
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
#  ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % ad_tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % ad_tflux(i,j,k,2,itrc) = IniVal
              END DO
#  endif
              FORCES(ng) % ad_stflx(i,j,itrc) = IniVal
              FORCES(ng) % ad_btflx(i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
      IF (model.eq.0) THEN
#  ifdef ADJUST_WSTRESS
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FORCES(ng) % b_sustr(i,j) = IniVal
            FORCES(ng) % b_svstr(i,j) = IniVal
#   ifdef FOUR_DVAR
            FORCES(ng) % e_sustr(i,j) = IniVal
            FORCES(ng) % e_svstr(i,j) = IniVal
            DO k=1,Nfrec(ng)
              FORCES(ng) % d_sustr(i,j,k) = IniVal
              FORCES(ng) % d_svstr(i,j,k) = IniVal
            END DO
#   endif
          END DO
        END DO
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
        DO itrc=1,NT(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FORCES(ng) % b_stflx(i,j,itrc) = IniVal
#   ifdef FOUR_DVAR
              FORCES(ng) % e_stflx(i,j,itrc) = IniVal
              DO k=1,Nfrec(ng)
                FORCES(ng) % d_stflx(i,j,k,itrc) = IniVal
              END DO
#   endif
            END DO
          END DO
        END DO
#  endif
      END IF
# endif
#endif
!
      RETURN
      END SUBROUTINE initialize_forces

      END MODULE mod_forces
