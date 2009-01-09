#include "cppdefs.h"
      MODULE mod_forces
!
!svn $Id: mod_forces.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Surface momentum stresses.                                          !
!                                                                      !
!  sustr        Kinematic surface momentum flux (wind stress) in       !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  sustrG       Latest two-time snapshots of input "sustr" grided      !
!                 data used for interpolation.                         !
!  svstr        Kinematic surface momentum flux (wind stress) in       !
!                 the ETA-direction (m2/s2) at horizontal V-points.    !
!  svstrG       Latest two-time snapshots of input "svstr" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom momentum stresses.                                           !
!                                                                      !
!  bustr        Kinematic bottom momentum flux (bottom stress) in      !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  bvstr        Kinematic bottom momentum flux (bottom stress) in      !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!                                                                      !
!  Surface wind induced waves.                                         !
!                                                                      !
!  Hwave        Surface wind induced wave height (m).                  !
!  HwaveG       Latest two-time snapshots of input "Hwave" grided      !
!                 data used for interpolation.                         !
!  Dwave        Surface wind induced wave direction (radians).         !
!  DwaveG       Latest two-time snapshots of input "Dwave" grided      !
!                 data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output       !
!  LwaveG       Latest two-time snapshots of input "Lwave" grided      !
!                 data used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).                  !
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided  !
!                 data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).                   !
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided  !
!                 data used for interpolation.                         !
!  Ub_swan      Bottom orbital velocity read in from swan output       !
!  Ub_swanG     Latest two-time snapshots of input "Ub_swan" grided    !
!                 data used for interpolation.                         !
!  wave_dissip  Wave dissipation                                       !
!  wave_dissipG Latest two-time snapshots of input "wave_dissip"       !
!                 gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.    !
!  Wave_breakG  Latest two-time snapshots of input "wave_break"        !
!                 gridded data used for interpolation.                 !
!                                                                      !
!  Solar shortwave radiation flux.                                     !
!                                                                      !
!  srflx        Kinematic surface shortwave solar radiation flux       !
!                 (Celsius m/s) at horizontal RHO-points               !
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
!  lhflx        Kinematic net sensible heat flux (degC m/s).           !
!  lrflx        Kinematic net longwave radiation (degC m/s).           !
!  shflx        Kinematic net sensible heat flux (degC m/s).           !
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
!  stflx        Kinematic surface flux of tracer type variables        !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                 horizontal RHO-points.                               !
!  stflxG       Latest two-time snapshots of input "stflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflx        Kinematic bottom flux of tracer type variables         !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                horizontal RHO-points.                                !
!  btflxG       Latest two-time snapshots of input "btflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Kinematic surface net heat flux sensitivity to SST,    !
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

        implicit none

        TYPE T_FORCES
!
!  Nonlinear model state.
!
          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
#if !defined ANA_SMFLUX && !defined BULK_FLUXES
          real(r8), pointer :: sustrG(:,:,:)
          real(r8), pointer :: svstrG(:,:,:)
#endif
#ifdef ADJUST_WSTRESS
          real(r8), pointer :: ustr(:,:,:,:)
          real(r8), pointer :: vstr(:,:,:,:)
#endif
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)

#ifdef WAVES_DIR
          real(r8), pointer :: Dwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: DwaveG(:,:,:)
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

#if defined BBL_MODEL || defined WAVES_OCEAN
          real(r8), pointer :: Ub_swan(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Ub_swanG(:,:,:)
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN
          real(r8), pointer :: wave_dissip(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: wave_dissipG(:,:,:)
# endif
#endif

#if defined SVENDSEN_ROLLER
          real(r8), pointer :: wave_break(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: wave_breakG(:,:,:)
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
          real(r8), pointer :: srflx(:,:)
#  ifndef ANA_SRFLUX
          real(r8), pointer :: srflxG(:,:,:)
#  endif
# endif

# ifdef CLOUDS
          real(r8), pointer :: cloud(:,:)
#  ifndef ANA_CLOUD
          real(r8), pointer :: cloudG(:,:,:)
#  endif
# endif
# ifdef BULK_FLUXES
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lrflx(:,:)
#  ifndef LONGWAVE
          real(r8), pointer :: lrflxG(:,:,:)
#  endif
          real(r8), pointer :: shflx(:,:)
# endif

# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
          real(r8), pointer :: Hair(:,:)
#  ifndef ANA_HUMIDITY
          real(r8), pointer :: HairG(:,:,:)
#  endif
          real(r8), pointer :: Tair(:,:)
#  ifndef ANA_TAIR
          real(r8), pointer :: TairG(:,:,:)
#  endif
# endif

# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          real(r8), pointer :: Pair(:,:)
#  ifndef ANA_PAIR
          real(r8), pointer :: PairG(:,:,:)
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM
          real(r8), pointer :: Uwind(:,:)
          real(r8), pointer :: Vwind(:,:)
#  ifndef ANA_WINDS
          real(r8), pointer :: UwindG(:,:,:)
          real(r8), pointer :: VwindG(:,:,:)
#  endif
# endif

# ifdef BULK_FLUXES
          real(r8), pointer :: rain(:,:)
#  ifndef ANA_RAIN
          real(r8), pointer :: rainG(:,:,:)
#  endif
#  ifdef EMINUSP
          real(r8), pointer :: evap(:,:)
#  endif
# endif

          real(r8), pointer :: stflx(:,:,:)
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
          real(r8), pointer :: stflxG(:,:,:,:)
# endif
# ifdef ADJUST_STFLUX
          real(r8), pointer :: tflux(:,:,:,:,:)
# endif
          real(r8), pointer :: btflx(:,:,:)
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
          real(r8), pointer :: btflxG(:,:,:,:)
# endif

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

# ifdef ECOSIM
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
          real(r8), pointer :: tl_stflx(:,:,:)
          real(r8), pointer :: tl_btflx(:,:,:)
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

        TYPE (T_FORCES), allocatable :: FORCES(:)

      CONTAINS

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
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
!  Nonlinear model state
!
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )

#if !defined ANA_SMFLUX && !defined BULK_FLUXES
      allocate ( FORCES(ng) % sustrG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % svstrG(LBi:UBi,LBj:UBj,2) )
#endif
#ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      allocate ( FORCES(ng) % vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
#endif
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )

#ifdef WAVES_DIR
      allocate ( FORCES(ng) % Dwave(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % DwaveG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef WAVES_HEIGHT
      allocate ( FORCES(ng) % Hwave(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % HwaveG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef WAVES_LENGTH
      allocate ( FORCES(ng) % Lwave(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % LwaveG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef WAVES_TOP_PERIOD
      allocate ( FORCES(ng) % Pwave_top(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_topG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef WAVES_BOT_PERIOD
      allocate ( FORCES(ng) % Pwave_bot(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_botG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN
      allocate ( FORCES(ng) % Ub_swan(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Ub_swanG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN
      allocate ( FORCES(ng) % wave_dissip(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_dissipG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#if defined SVENDSEN_ROLLER
      allocate ( FORCES(ng) % wave_break(LBi:UBi,LBj:UBj) )
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_breakG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SRFLUX
      allocate ( FORCES(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef CLOUDS
      allocate ( FORCES(ng) % cloud(LBi:UBi,LBj:UBj) )

#  ifndef ANA_CLOUD
      allocate ( FORCES(ng) % cloudG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef BULK_FLUXES
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )

#  ifndef LONGWAVE
      allocate ( FORCES(ng) % lrflxG(LBi:UBi,LBj:UBj,2) )
#  endif

      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
# endif

# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      allocate ( FORCES(ng) % Hair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_HUMIDITY
      allocate ( FORCES(ng) % HairG(LBi:UBi,LBj:UBj,2) )
#  endif

      allocate ( FORCES(ng) % Tair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_TAIR
      allocate ( FORCES(ng) % TairG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      allocate ( FORCES(ng) % Pair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_PAIR
      allocate ( FORCES(ng) % PairG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined BULK_FLUXES || defined ECOSIM
      allocate ( FORCES(ng) % Uwind(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Vwind(LBi:UBi,LBj:UBj) )

#  ifndef ANA_WINDS
      allocate ( FORCES(ng) % UwindG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % VwindG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef BULK_FLUXES
      allocate ( FORCES(ng) % rain(LBi:UBi,LBj:UBj) )

#  ifndef ANA_RAIN
      allocate ( FORCES(ng) % rainG(LBi:UBi,LBj:UBj,2) )
#  endif
#  ifdef EMINUSP
      allocate ( FORCES(ng) % evap(LBi:UBi,LBj:UBj) )
#  endif
# endif

      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )

# if !defined ANA_STFLUX || !defined ANA_SSFLUX
      allocate ( FORCES(ng) % stflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
# endif
# ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tflux(LBi:UBi,LBj:UBj,nfrec(ng),          &
     &                              2,NT(ng)) )
# endif
      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )

# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
      allocate ( FORCES(ng) % btflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
# endif

# ifdef QCORRECTION
      allocate ( FORCES(ng) % dqdt(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % sst(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SST
      allocate ( FORCES(ng) % dqdtG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % sstG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      allocate ( FORCES(ng) % sss(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SSS
      allocate ( FORCES(ng) % sssG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef ECOSIM
      allocate ( FORCES(ng) % SpecIr(LBi:UBi,LBj:UBj,NBands) )
      allocate ( FORCES(ng) % avcos(LBi:UBi,LBj:UBj,NBands) )
# endif

#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state
!
      allocate ( FORCES(ng) % tl_sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % tl_svstr(LBi:UBi,LBj:UBj) )
# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      allocate ( FORCES(ng) % tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
# endif
      allocate ( FORCES(ng) % tl_bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % tl_bvstr(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( FORCES(ng) % tl_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( FORCES(ng) % tl_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tl_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % tl_srflx(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % tl_lhflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % tl_lrflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % tl_shflx(LBi:UBi,LBj:UBj) )
#   ifdef EMINUSP
      allocate ( FORCES(ng) % tl_evap(LBi:UBi,LBj:UBj) )
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state
!
      allocate ( FORCES(ng) % ad_sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_svstr(LBi:UBi,LBj:UBj) )
# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      allocate ( FORCES(ng) % ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
# endif
      allocate ( FORCES(ng) % ad_bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_bvstr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_bustr_sol(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_bvstr_sol(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( FORCES(ng) % ad_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( FORCES(ng) % ad_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % ad_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % ad_srflx(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % ad_lhflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_lrflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % ad_shflx(LBi:UBi,LBj:UBj) )
#   ifdef EMINUSP
      allocate ( FORCES(ng) % ad_evap(LBi:UBi,LBj:UBj) )
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
      allocate ( FORCES(ng) % b_svstr(LBi:UBi,LBj:UBj) )
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % b_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
#  endif
# endif
# ifdef FOUR_DVAR
#  ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      allocate ( FORCES(ng) % d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      allocate ( FORCES(ng) % e_sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % e_svstr(LBi:UBi,LBj:UBj) )
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % d_stflx(LBi:UBi,LBj:UBj,                  &
     &                                Nfrec(ng),NT(ng)) )
      allocate ( FORCES(ng) % e_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE allocate_forces

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
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
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
#if !defined ANA_SMFLUX && !defined BULK_FLUXES
            FORCES(ng) % sustrG(i,j,1) = IniVal
            FORCES(ng) % sustrG(i,j,2) = IniVal
            FORCES(ng) % svstrG(i,j,1) = IniVal
            FORCES(ng) % svstrG(i,j,2) = IniVal
#endif
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal
#ifdef WAVES_DIR
            FORCES(ng) % Dwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % DwaveG(i,j,1) = IniVal
            FORCES(ng) % DwaveG(i,j,2) = IniVal
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
            FORCES(ng) % Ub_swan(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Ub_swanG(i,j,1) = IniVal
            FORCES(ng) % Ub_swanG(i,j,2) = IniVal
# endif
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN
            FORCES(ng) % Wave_dissip(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_dissipG(i,j,1) = IniVal
            FORCES(ng) % Wave_dissipG(i,j,2) = IniVal
# endif
#endif
#if defined SVENDSEN_ROLLER
            FORCES(ng) % Wave_break(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_breakG(i,j,1) = IniVal
            FORCES(ng) % Wave_breakG(i,j,2) = IniVal
# endif
#endif

#ifdef SOLVE3D
# ifdef SHORTWAVE
            FORCES(ng) % srflx(i,j) = IniVal
#  ifndef ANA_SRFLUX
            FORCES(ng) % srflxG(i,j,1) = IniVal
            FORCES(ng) % srflxG(i,j,2) = IniVal
#  endif
# endif
# ifdef CLOUDS
            FORCES(ng) % cloud(i,j) = IniVal
#  ifndef ANA_CLOUD
            FORCES(ng) % cloudG(i,j,1) = IniVal
            FORCES(ng) % cloudG(i,j,2) = IniVal
#  endif
# endif
# ifdef BULK_FLUXES
            FORCES(ng) % lhflx(i,j) = IniVal
            FORCES(ng) % lrflx(i,j) = IniVal
            FORCES(ng) % shflx(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
            FORCES(ng) % Hair(i,j) = IniVal
            FORCES(ng) % Tair(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            FORCES(ng) % Pair(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM
            FORCES(ng) % Uwind(i,j) = IniVal
            FORCES(ng) % Vwind(i,j) = IniVal
# endif
# ifdef BULK_FLUXES
            FORCES(ng) % rain(i,j) = IniVal
#  ifdef EMINUSP
            FORCES(ng) % evap(i,j) = IniVal
#  endif
# endif
# if !defined LONGWAVE && defined BULK_FLUXES
            FORCES(ng) % lrflxG(i,j,1) = IniVal
            FORCES(ng) % lrflxG(i,j,2) = IniVal
# endif
# if defined BULK_FLUXES || defined ECOSIM
#  ifndef ANA_HUMIDITY
            FORCES(ng) % HairG(i,j,1) = IniVal
            FORCES(ng) % HairG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
#  ifndef ANA_PAIR
            FORCES(ng) % PairG(i,j,1) = IniVal
            FORCES(ng) % PairG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM
#  ifndef ANA_TAIR
            FORCES(ng) % TairG(i,j,1) = IniVal
            FORCES(ng) % TairG(i,j,2) = IniVal
#  endif
#  ifndef ANA_WINDS
            FORCES(ng) % UwindG(i,j,1) = IniVal
            FORCES(ng) % UwindG(i,j,2) = IniVal
            FORCES(ng) % VwindG(i,j,1) = IniVal
            FORCES(ng) % VwindG(i,j,2) = IniVal
#  endif
# endif
# if !defined ANA_RAIN && defined BULK_FLUXES
            FORCES(ng) % rainG(i,j,1) = IniVal
            FORCES(ng) % rainG(i,j,2) = IniVal
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
# ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % tflux(i,j,k,2,itrc) = IniVal
              END DO
# endif
              FORCES(ng) % stflx(i,j,itrc) = IniVal
              FORCES(ng) % btflx(i,j,itrc) = IniVal
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
              FORCES(ng) % stflxG(i,j,1,itrc) = IniVal
              FORCES(ng) % stflxG(i,j,2,itrc) = IniVal
# endif
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
              FORCES(ng) % btflxG(i,j,1,itrc) = IniVal
              FORCES(ng) % btflxG(i,j,2,itrc) = IniVal
# endif
            END DO
# ifdef ECOSIM
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
              FORCES(ng) % tl_stflx(i,j,itrc) = IniVal
              FORCES(ng) % tl_btflx(i,j,itrc) = IniVal
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

      RETURN
      END SUBROUTINE initialize_forces

      END MODULE mod_forces
