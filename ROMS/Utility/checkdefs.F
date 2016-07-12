#include "cppdefs.h"
      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ibbl = 0
      integer :: ibiology = 0
      integer :: idriver = 0
      integer :: itrcHadv = 0
      integer :: itrcVadv = 0
      integer :: itrcHadvtl = 0
      integer :: itrcVadvtl = 0
      integer :: ivelHadv = 0
      integer :: ivelVadv = 0
      integer :: ivmix = 0
      integer :: nearshore = 0
#ifdef SEDIMENT
      integer :: sbed_load = 0
      integer :: sbed_type = 0
#endif
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t22,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
#if defined AD_AVERAGES && defined ADJOINT
      IF (Master) WRITE (stdout,20) 'AD_AVERAGES',                      &
     &   'Writing out time-averaged adjoint model fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' AD_AVERAGES,'
#endif
#if defined ADD_FSOBC && defined SSH_TIDES
      IF (Master) WRITE (stdout,20) 'ADD_FSOBC',                        &
     &   'Adding tidal elevation to processed OBC data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_FSOBC,'
#endif
#if defined ADD_M2OBC && defined UV_TIDES
      IF (Master) WRITE (stdout,20) 'ADD_M2OBC',                        &
     &   'Adding tidal currents to processed OBC data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_M2OBC,'
#endif
#ifdef ADJOINT
      IF (Master) WRITE (stdout,20) 'ADJOINT',                          &
     &   'Adjoint Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ADJOINT,'
#endif
#if defined ADJUST_BOUNDARY     && \
   (defined CORRELATION         || defined IS4DVAR           || \
    defined IS4DVAR_SENSITIVITY || defined SENSITIVITY_4DVAR || \
    defined W4DPSAS             || defined W4DVAR)
      IF (Master) WRITE (stdout,20) 'ADJUST_BOUNDARY',                  &
     &   'Including boundary conditions in 4DVar state estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ADJUST_BOUNDARY,'
#endif
#if defined ADJUST_STFLUX       && \
   (defined CORRELATION         || defined IS4DVAR           || \
    defined IS4DVAR_SENSITIVITY || defined SENSITIVITY_4DVAR || \
    defined W4DPSAS             || defined W4DVAR)
      IF (Master) WRITE (stdout,20) 'ADJUST_STFLUX',                    &
     &   'Including surface tracer flux in 4DVar state estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ADJUST_STFLUX,'
#endif
#if defined ADJUST_WSTRESS      && \
   (defined CORRELATION         || defined IS4DVAR           || \
    defined IS4DVAR_SENSITIVITY || defined SENSITIVITY_4DVAR || \
    defined W4DPSAS             || defined W4DVAR)
      IF (Master) WRITE (stdout,20) 'ADJUST_WSTRESS',                   &
     &   'Including surface wind stress in 4DVar state estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ADJUST_WSTRESS,'
#endif
#if defined AD_IMPULSE
      IF (Master) WRITE (stdout,20) 'AD_IMPULSE',                       &
     &   'Force adjoint model with intermittent impulses.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' AD_IMPULSE,'
#endif
#ifdef ADM_DRIVER
      IF (Master) WRITE (stdout,20) 'ADM_DRIVER',                       &
     &   'Generic adjoint model driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ADM_DRIVER,'
      idriver=idriver+1
#endif
#ifdef AD_SENSITIVITY
      IF (Master) WRITE (stdout,20) 'AD_SENSITIVITY',                   &
     &   'Adjoint Sensitivity Analysis.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' AD_SENSITIVITY,'
      idriver=idriver+1
#endif
#ifdef AFT_EIGENMODES
      IF (Master) WRITE (stdout,20) 'AFT_EIGENMODES',                   &
     &   'Adjoint Finite Time Eigenvalues.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' AFT_EIGENMODES,'
      idriver=idriver+1
#endif
#if defined AGE_MEAN && defined T_PASSIVE && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'AGE_MEAN',                         &
     &   'Computing Mean Age of inert passive tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AGE_MEAN,'
      idriver=idriver+1
#endif
#ifdef AICLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'AICLIMATOLOGY',                    &
     &   'Ice concentration climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' AICLIMATOLOGY,'
#endif
#ifdef AICLM_NUDGING
      IF (Master) WRITE (stdout,20) 'AICLM_NUDGING',                    &
     &   'Nudging to ice concentration climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' AICLM_NUDGING,'
#endif
#if defined ALBEDO && defined ANA_SRFLUX
      IF (Master) WRITE (stdout,20) 'ALBEDO',                           &
     &   'Shortwave radiation with albedo correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ALBEDO,'
#endif
#if defined ALBEDO_CLOUD
      IF (Master) WRITE (stdout,20) 'ALBEDO_CLOUD',                     &
     &   'Albedo as a function of cloudiness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ALBEDO_CLOUD,'
#endif
#if defined ALBEDO_CSIM
      IF (Master) WRITE (stdout,20) 'ALBEDO_CSIM',                      &
     &   'Albedo computed according to CSIM3.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ALBEDO_CSIM,'
#endif
#if defined ALBEDO_CURVE
      IF (Master) WRITE (stdout,20) 'ALBEDO_CURVE',                     &
     &   'Shortwave radiation adjusted by albedo curve.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ALBEDO_CURVE,'
#endif
#if defined ALBEDO_FILE
      IF (Master) WRITE (stdout,20) 'ALBEDO_FILE',                      &
     &   'Albedo read from a forcing file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ALBEDO_FILE,'
#endif
#ifdef ANA_ALBEDO
      IF (Master) WRITE (stdout,20) 'ANA_ALBEDO',                       &
     &   'Analytical surface albedo.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_ALBEDO,'
#endif
#if defined BIOLOGY && defined ANA_BIOLOGY
      IF (Master) WRITE (stdout,20) 'ANA_BIOLOGY',                      &
     &   'Analytical biology initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_BIOLOGY,'
#endif
#if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
# ifdef ANA_BPFLUX
      IF (Master) WRITE (stdout,20) 'ANA_BPFLUX',                       &
     &   'Analytical bottom passive tracers fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BPFLUX,'
# endif
#endif
#ifdef ANA_BSFLUX
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
#endif
#ifdef ANA_BTFLUX
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
#endif
#if defined ANA_CLOUD && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ANA_CLOUD',                        &
     &   'Analytical cloud fraction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_CLOUD,'
#endif
#if defined ANA_DRAG && defined UV_DRAG_GRID
      IF (Master) WRITE (stdout,20) 'ANA_DRAG_GRID',                    &
# if defined UV_LOGDRAG
     &   'Analytical spatially varying bottom roughness length.'
# elif defined UV_LDRAG
     &   'Analytical spatially varying linear drag coefficient.'
# elif defined UV_QDRAG
     &   'Analytical spatially varying quadratic drag coefficient.'
# endif
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_DRAG,'
#endif
#ifdef ANA_DIAG
      IF (Master) WRITE (stdout,20) 'ANA_DIAG',                         &
     &   'Customized diagnostics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_DIAG,'
#endif
#if defined ANA_DQDSST && defined QCORRECTION
      IF (Master) WRITE (stdout,20) 'ANA_DQDSST',                          &
     &   'Analytical surface heat flux sensitivity to SST.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_DQDSST,'
#endif
#ifdef ANA_FSOBC
      IF (Master) WRITE (stdout,20) 'ANA_FSOBC',                        &
     &   'Analytical free-surface boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_FSOBC,'
#endif
#ifdef ANA_GRID
      IF (Master) WRITE (stdout,20) 'ANA_GRID',                         &
     &   'Analytical grid set-up.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_GRID,'
#endif
#if defined ANA_HUMIDITY && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ANA_HUMIDITY',                     &
     &   'Analytical surface air humidity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_HUMIDITY,'
#endif
#ifdef ANA_ICE
      IF (Master) WRITE (stdout,20) 'ANA_ICE',                          &
     &   'Analytical ice initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_ICE,'
#endif
#ifdef ANA_INITIAL
      IF (Master) WRITE (stdout,20) 'ANA_INITIAL',                      &
     &   'Analytical initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_INITIAL,'
#endif
#ifdef ANA_LRFLUX
      IF (Master) WRITE (stdout,20) 'ANA_LRFLUX',                       &
     &   'Analytical incoming longwave radiation flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_LRFLUX,'
#endif
#ifdef ANA_M2CLIMA
      IF (Master) WRITE (stdout,20) 'ANA_M2CLIMA',                      &
     &   'Analytical 2D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_M2CLIMA,'
#endif
#ifdef ANA_M2OBC
      IF (Master) WRITE (stdout,20) 'ANA_M2OBC',                        &
     &   'Analytical 2D momentum boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_M2OBC,'
#endif
#if defined ANA_M3CLIMA && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_M3CLIMA',                      &
     &   'Analytical 3D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_M3CLIMA,'
#endif
#if defined ANA_M3OBC && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_M3OBC',                        &
     &   'Analytical 3D momentum boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_M3OBC,'
#endif
#ifdef ANA_MASK
      IF (Master) WRITE (stdout,20) 'ANA_MASK',                         &
     &   'Analytical Land/Sea Masking.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_MASK,'
#endif
#ifdef ANA_NUDGCOEF
      IF (ANY(Lnudging(:))) THEN
        IF (Master) WRITE (stdout,20) 'ANA_NUDGCOEF',                   &
     &     'Analytical spatially varying nudging time-scales.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+14)=' ANA_NUDGCOEF,'
      END IF
#endif
#if defined ANA_PAIR && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ANA_PAIR',                         &
     &   'Analytical surface air pressure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_PAIR,'
#endif
#if defined ANA_PASSIVE && defined T_PASSIVE && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_PASSIVE',                      &
     &   'Analytical initial conditions for inert tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_PASSIVE,'
#endif
#if defined ANA_PERTURB
      IF (Master) WRITE (stdout,20) 'ANA_PERTURB',                      &
     &   'Perturb initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_PERTURB,'
#endif
#ifdef ANA_PSOURCE
      IF (Master) WRITE (stdout,20) 'ANA_PSOURCE',                      &
     &   'Analytical point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_PSOURCE,'
#endif
#if defined ANA_PTOBC && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_PTOBC',                        &
     &   'Analytical passive tracers boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_PTOBC,'
#endif
#if defined ANA_RAIN && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ANA_RAIN',                         &
     &   'Analytical rain fall rate.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_RAIN,'
#endif
#if defined ANA_SCOPE        && \
   (defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI)
      IF (Master) WRITE (stdout,20) 'ANA_SCOPE',                        &
     &   'Analytical spatial scope masking arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SCOPE,'
#endif
#if defined ANA_SEDIMENT && (defined SEDIMENT || defined BBL_MODEL)
      IF (Master) WRITE (stdout,20) 'ANA_SEDIMENT',                     &
     &   'Analytical sediment initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_SEDIMENT,'
#endif
#if defined VEGETATION && defined ANA_VEGETATION 
      IF (Master) WRITE (stdout,20) 'ANA_VEGETATION',                   &
     &   'Analytical vegetation initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ANA_VEGETATION,'
#endif
#ifdef ANA_SMFLUX
      IF (Master) WRITE (stdout,20) 'ANA_SMFLUX',                       &
     &   'Analytical kinematic surface momentum flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SMFLUX,'
#endif
#if defined ANA_SNOW
      IF (Master) WRITE (stdout,20) 'ANA_SNOW',                         &
     &   'Analytical snowfall rate.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_SNOW,'
#endif
#if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
# ifdef ANA_SPFLUX
      IF (Master) WRITE (stdout,20) 'ANA_SPFLUX',                       &
     &   'Analytical surface passive tracer fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SPFLUX,'
# endif
#endif
#ifdef ANA_SPINNING
      IF (Master) WRITE (stdout,20) 'ANA_SPINNING',                     &
     &   'Analytical time-varying rotation force.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_SPINNING,'
#endif
#ifdef ANA_SPONGE
      IF (ANY(Lsponge)) THEN
        IF (Master) WRITE (stdout,20) 'ANA_SPONGE',                     &
     &     'Analytical enhanced viscosity/diffusion sponge.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+12)=' ANA_SPONGE,'
      END IF
#endif
#ifdef ANA_SRFLUX
      IF (Master) WRITE (stdout,20) 'ANA_SRFLUX',                       &
     &   'Analytical kinematic shortwave radiation flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SRFLUX,'
#endif
#ifdef ANA_SSH
      IF (Master) WRITE (stdout,20) 'ANA_SSH',                          &
     &   'Analytical sea surface height climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SSH,'
#endif
#ifdef ANA_SSS
      IF (Master) WRITE (stdout,20) 'ANA_SSS',                          &
     &   'Analytical sea surface salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SSS,'
#endif
#if defined ANA_SST && defined QCORRECTION
      IF (Master) WRITE (stdout,20) 'ANA_SST',                          &
     &   'Analytical sea surface temperature, SST.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SST,'
#endif
#ifdef ANA_SSFLUX
      IF (Master) WRITE (stdout,20) 'ANA_SSFLUX',                       &
     &   'Analytical kinematic surface salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SSFLUX,'
#endif
#ifdef ANA_STFLUX
      IF (Master) WRITE (stdout,20) 'ANA_STFLUX',                       &
     &   'Analytical kinematic surface temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_STFLUX,'
#endif
#if defined ANA_TCLIMA && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_TCLIMA',                       &
     &   'Analytical tracer climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_TCLIMA,'
#endif
#if defined ANA_TOBC && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_TOBC',                         &
     &   'Analytical tracers boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_TOBC,'
#endif
#if defined ANA_VMIX && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'ANA_VMIX',                         &
     &   'Analytical vertical mixing coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_VMIX,'
      ivmix=ivmix+1
#endif
#if defined ANA_WINDS && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ANA_WINDS',                        &
     &   'Analytical surface wind components.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WINDS,'
#endif
#if defined ANA_WTYPE && defined WTYPE_GRID    && \
   (defined LMD_SKPP  || defined SOLAR_SOURCE)
      IF (Master) WRITE (stdout,20) 'ANA_WTYPE',                        &
     &   'Analytical spatially varying Jerlov water type index.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WTYPE,'
#endif
#ifdef ANA_WWAVE
      IF (Master) WRITE (stdout,20) 'ANA_WWAVE',                        &
     &   'Analytical wind induced waves.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WWAVE,'
#endif
#ifdef ARRAY_MODES
      IF (Master) WRITE (stdout,20) 'ARRAY_MODES',                      &
     &   'Stabilized Representer Matrix Array Modes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ARRAY_MODES,'
      idriver=idriver+1
# ifdef ARRAY_MODES_SPLIT
      IF (Master) WRITE (stdout,20) 'ARRAY_MODES_SPLIT',                &
     &   'Separate analysis due to IC, surface forcing, and OBC'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' ARRAY_MODES_SPLIT,'
# endif
#endif
#ifdef ASSUMED_SHAPE
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
#endif
#ifdef ATM_PRESS
      IF (Master) WRITE (stdout,20) 'ATM_PRESS',                        &
     &   'Impose atmospheric pressure onto sea surface.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ATM_PRESS,'
#endif
#ifdef ATM2OCN_FLUXES
      IF (Master) WRITE (stdout,20) 'ATM2OCN_FLUXES',                   &
     &   'Surface fluxes from ATM Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ATM2OCN_FLUXES,'
#endif
#ifdef AVERAGES
      IF (Master) WRITE (stdout,20) 'AVERAGES',                         &
     &   'Writing out time-averaged nonlinear model fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AVERAGES,'
# ifdef AVERAGES2
      IF (Master) WRITE (stdout,20) 'AVERAGES2',                        &
     &   'Writing out secondary time-averaged fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' AVERAGES2,'
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      IF (Master) WRITE (stdout,20) 'AVERAGES_DETIDE',                  &
     &   'Writing out time-averaged nonlinear model detided fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' AVERAGES_DETIDE,'
# endif
# ifdef AVERAGES_UV_KIRBY
      IF (Master) WRITE (stdout,20) 'AVERAGES_UV_KIRBY',                &
     &   'Writing out time-averaged uwave and vwave Kirby velocities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' AVERAGES_UV_KIRBY,'
# endif
# ifdef AVERAGES_WEC
      IF (Master) WRITE (stdout,20) 'AVERAGES_WEC',                     &
     &   'Writing out time-averaged nearshore radiation fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' AVERAGES_WEC,'
# endif
#endif
#if defined BACKGROUND && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'BACKGROUND',                       &
     &   'Include background cost function.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BACKGROUND,'
#endif
#ifdef BALANCE_OPERATOR
      IF (Master) WRITE (stdout,20) 'BALANCE_OPERATOR',                 &
     &   'Error Covariance Multivariate Balance Operator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' BALANCE_OPERATOR,'
#endif
#if defined SEDIMENT && defined BEDLOAD_MPM
      IF (Master) WRITE (stdout,20) 'BEDLOAD_MPM',                      &
     &   'Activate bed load sediment transport Meyer-Peter-Mueller.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BEDLOAD_MPM,'
      sbed_load=sbed_load+1
#endif
#if defined SEDIMENT && defined BEDLOAD_SOULSBY
      IF (Master) WRITE (stdout,20) 'BEDLOAD_SOULSBY',                  &
     &   'Activate bed load sediment transport Soulsby formula.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' BEDLOAD_SOULSBY,'
#endif
#if defined BGQC && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'BGQC',                             &
     &   'Background quality control of observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' BGQC,'
#endif
#if defined BEOFS_ONLY && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'BEOFS_ONLY',                       &
     &   'Compute EOFs of the background error covariance matrix.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BEOFS_ONLY,'
#endif
#ifdef BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'BIO_FENNEL',                       &
     &   'Fennel et al. (2006) nitrogen-based model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BIO_FENNEL,'
      ibiology=ibiology+1
#endif
#if defined BNORM && defined HESSIAN_SV
      IF (Master) WRITE (stdout,20) 'BNORM',                            &
     &   'Background normalization of Hessian singular vectors.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+7)=' BNORM,'
#endif
#ifdef ESTUARYBGC 
      IF (Master) WRITE (stdout,20) 'ESTUARYBGC',                       &
     &   'Estuary bio-geo-chemical model (Testa et. al).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ESTUARYBGC,'
      ibiology=ibiology+1
#endif
#if defined BIO_SEDIMENT && \
   (defined BIO_FENNEL   || defined NEMURO || defined ECOSIM)
      IF (Master) WRITE (stdout,20) 'BIO_SEDIMENT',                     &
     &   'Restore fallen particulate material to the nutrient pool.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' BIO_SEDIMENT,'
#endif
#ifdef BIO_UMAINE
      IF (Master) WRITE (stdout,20) 'BIO_UMAINE',                       &
     &   'UMaine CoSiNE ecosystem model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BIO_UMAINE,'
      ibiology=ibiology+1
#endif
#ifdef BODYFORCE
      IF (Master) WRITE (stdout,20) 'BODYFORCE',                        &
     &   'Momentum stresses as body-forces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' BODYFORCE,'
#endif
#ifdef BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES',                      &
     &   'Surface bulk fluxes parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
#endif
#ifdef BULK_FLUXES2D
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES2D',                    &
     &   'Surface bulk fluxes parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' BULK_FLUXES2D,'
#endif
#ifdef BVF_MIXING
      IF (Master) WRITE (stdout,20) 'BVF_MIXING',                       &
     &   'Brunt-Vaisala frequency based vertical mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BVF_MIXING,'
      ivmix=ivmix+1
#endif
#if defined CANUTO_A && defined GLS_MIXING
      IF (Master) WRITE (stdout,20) 'CANUTO_A',                         &
     &   'Canuto A-stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CANUTO_A,'
#endif
#if defined CANUTO_B && defined GLS_MIXING
      IF (Master) WRITE (stdout,20) 'CANUTO_B',                         &
     &   'Canuto B-stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CANUTO_B,'
#endif
#if defined CARBON && defined BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'CARBON',                           &
     &   'Add Carbon constituents.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CARBON,'
#endif
#if defined CARBON && defined ESTUARYBGC
      IF (Master) WRITE (stdout,20) 'CARBON',                           &
     &   'Add Carbon constituents.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CARBON,'
#endif
#ifdef CCSM_FLUXES
      IF (Master) WRITE (stdout,20) 'CCSM_FLUXES',                      &
     &   'Surface bulk fluxes parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' CCSM_FLUXES,'
#endif
#ifdef CCSM_FLUXES2D
      IF (Master) WRITE (stdout,20) 'CCSM_FLUXES2D',                    &
     &   'Surface bulk fluxes parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' CCSM_FLUXES2D,'
#endif
#if defined CELERITY_READ && defined FORWARD_READ
      IF (Master) WRITE (stdout,20) 'CELERITY_READ',                    &
     &   'Read in radiation celerity from forward file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' CELERETY_READ,'
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
      IF (Master) WRITE (stdout,20) 'CELERITY_WRITE',                   &
     &   'Write out radiation celerity in forward file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' CELERITY_WRITE,'
#endif
#ifdef CICE_MODEL
      IF (Master) WRITE (stdout,20) 'CICE_MODEL',                       &
     &   'Include CICE Ice Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' CICE_MODEL,'
#endif
#ifdef CLIPPING
      IF (Master) WRITE (stdout,20) 'CLIPPING',                         &
     &   'Stabilized Representer Matrix clipping spectrum analysis.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CLIPPING,'
      idriver=idriver+1
# ifdef CLIPPING_SPLIT
      IF (Master) WRITE (stdout,20) 'CLIPPING_SPLIT',                   &
     &   'Separate analysis due to IC, surface forcing, and OBC'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' CLIPPING_SPLIT,'
# endif
#endif
#ifdef CORRELATION
      IF (Master) WRITE (stdout,20) 'CORRELATION',                      &
     &   'Background-error Correlation Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' CORRELATION,'
      idriver=idriver+1
#endif
#if defined DATALESS_LOOPS && defined W4DVAR
      IF (Master) WRITE (stdout,20) 'DATALESS_LOOPS',                   &
     &   'Testing convergence of RPM Picard iterations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DATALESS_LOOPS,'
#endif
#if defined GLS_MIXING && defined CHARNOK
      IF (Master) WRITE (stdout,20) 'CHARNOK',                          &
     &   'Charnok surface roughness from wind stress.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' CHARNOK,'
#endif
#if defined PROPAGATOR && defined CHECKPOINTING
      IF (Master) WRITE (stdout,20) 'CHECKPOINTING',                    &
     &   'Processing checkpointing NetCDF for GST analysis.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' CHECKPOINTING,'
#endif
#if defined CLIMA_TS_MIX && defined SOLVE3D && \
   (defined TS_DIF2      || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'CLIMA_TS_MIX',                     &
     &   'Horizontal diffusion of tracer perturbation (T-Tclm).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CLIMA_TS_MIX,'
#endif
#if defined GLS_MIXING && defined CRAIG_BANNER
      IF (Master) WRITE (stdout,20) 'CRAIG_BANNER',                     &
     &   'Craig and Banner wave breaking surface flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CRAIG_BANNER,'
#endif
#if defined COARE_OOST
      IF (Master) WRITE (stdout,20) 'COARE_OOST',                       &
     &   'Oost et al (2002) relation for ZoW in bulk fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' COARE_OOST,'
#endif
#if defined COARE_TAYLOR_YELLAND
      IF (Master) WRITE (stdout,20) 'COARE_TAYLOR_YELLAND',             &
     &   'Taylor and Yelland (2001) relation for ZoW in bulk fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+22)=' COARE_TAYLOR_YELLAND,'
#endif
#if defined COOL_SKIN && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'COOL_SKIN',                        &
     &   'Surface cool skin correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' COOL_SKIN,'
#endif
#ifdef CORE_FORCING
      IF (Master) WRITE (stdout,20) 'CORE_FORCING',                      &
     &   'Humidity from CORE forcing files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CORE_FORCING,'
#endif
#if defined COSINE2 && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'COSINE2',                          &
     &   'Cosine-squared shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' COSINE2,'
#endif
#ifdef CURVGRID
      IF (Master) WRITE (stdout,20) 'CURVGRID',                         &
     &   'Orthogonal curvilinear grid.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CURVGRID,'
#endif
#if defined DAILY_SHORTWAVE && defined RED_TIDE
      IF (Master) WRITE (stdout,20) 'DAILY_SHORTWAVE',                  &
     &   'Daily-averaged surface solar shortwave radiation flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DAILY_SHORTWAVE,'
#endif
#ifdef DEBUGGING
      IF (Master) WRITE (stdout,20) 'DEBUGGING',                        &
     &   'Internal debugging switch activated.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DEBUGGING,'
#endif
#ifdef DEEPWATER_WAVES
      IF (Master) WRITE (stdout,20) 'DEEPWATER_WAVES',                  &
     &   'Deep water waves approx in bulk fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DEEPWATER_WAVES,'
#endif
#if defined DEFLATE && defined HDF5
      IF (Master) WRITE (stdout,20) 'DEFLATE',                          &
     &   'Setting compression in output NetCDF-4/HDF5 files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DEFLATE,'
#endif
#if defined DENITRIFICATION && defined BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'DENITRIFICATION',                  &
     &   'Add denitrification processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DENITRIFICATION,'
#endif
#if defined DENITRIFICATION && defined ESTUARYBGC 
      IF (Master) WRITE (stdout,20) 'DENITRIFICATION',                  &
     &   'Add denitrification processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DENITRIFICATION,'
#endif
#if defined DIAGNOSTICS_BIO && defined BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_BIO',                  &
     &   'Computing and writing biological diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DIAGNOSTICS_BIO,'
#endif
#if defined DIAGNOSTICS_BIO && defined ESTUARYBGC
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_BIO',                  &
     &   'Computing and writing biological diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DIAGNOSTICS_BIO,'
#endif
#if defined SPECTRAL_LIGHT && defined BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'SPECTRAL_LIGHT',                   &
     &   'Computing and writing spectral light terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SPECTRAL_LIGHT,'
#endif
#if defined SPECTRAL_LIGHT && defined ESTUARYBGC
      IF (Master) WRITE (stdout,20) 'SPECTRAL_LIGHT',                   &
     &   'Computing and writing spectral light terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SPECTRAL_LIGHT,'
#endif
#if defined SPECTRAL_LIGHT && defined CDOM_DEFAULT
      IF (Master) WRITE (stdout,20) 'CDOM_DEFAULT',                     &
     &   'Constant default CDOM value for spectral light terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CDOM_DEFAULT,'
#endif
#if defined SPECTRAL_LIGHT && defined CDOM_VARIABLE
      IF (Master) WRITE (stdout,20) 'CDOM_VARIABLE',                    &
     &   'Variable CDOM for spectral light terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' CDOM_VARIABLE,'
#endif
#if defined SPECTRAL_LIGHT && defined MOD_SWR_SPECT
      IF (Master) WRITE (stdout,20) 'MOD_SWR_SPECT',                    &
     &   'modulate shortwave radiation PARsur using Gallegos spectrum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' MOD_SWR_SPECT,'
#endif
#if defined SPECTRAL_LIGHT && defined MOD_SWR_HOMO
      IF (Master) WRITE (stdout,20) 'MOD_SWR_HOMO',                     &
     &   'modulate shortwave radiation PARsur evenly through spectrum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' MOD_SWR_HOMO,'
#endif
#if defined SPECTRAL_LIGHT && defined CHL_BACKSCAT
      IF (Master) WRITE (stdout,20) 'CHL_BACKSCAT',                     &
     &   'Chlorophyll backscatter from ecosim (not in gallegos).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CHL_BACKSCAT,'
#endif
#if defined VEGETATION 
      IF (Master) WRITE (stdout,20) 'VEGETATION',                       &
     &   'Active submerged/emergent vegetation effects.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VEGETATION,'
#endif
#if defined VEGETATION && defined MARSH_WAVE_THRUST 
      IF (Master) WRITE (stdout,20) 'MARSH_WAVE_THRUST',                &
     &   'Vegetation module containing wave thrust effect on marshes .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' MARSH_WAVE_THRUST,'
#endif
#ifdef DIAGNOSTICS_TS
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_TS',                   &
     &   'Computing and writing tracer diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DIAGNOSTICS_TS,'
#endif
#ifdef DIAGNOSTICS_UV
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_UV',                   &
     &   'Computing and writing momentum diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DIAGNOSTICS_UV,'
#endif
#ifdef DIAPAUSE
      IF (Master) WRITE (stdout,20) 'DIAPAUSE',                   &
     &   'Floats moving vertically to simulate diapause.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' DIAPAUSE,'
#endif
#if defined DIFF_3DCOEF
      IF (Master) WRITE (stdout,20) 'DIFF_3DCOEF',                      &
     &   'Horizontal, time-dependent 3D diffusion coefficient.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' DIFF_3DCOEF,'
#endif
#if defined TS_DIF2 || defined TS_DIF4
# ifdef DIFF_GRID
      IF (Master) WRITE (stdout,20) 'DIFF_GRID',                        &
     &   'Horizontal diffusion coefficient scaled by grid size.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DIFF_GRID,'
# endif
#endif
#ifdef DIURNAL_SRFLUX
      IF (Master) WRITE (stdout,20) 'DIURNAL_SRFLUX',                   &
     &   'Modulate shortwave radiation by the local diurnal cycle.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DIURNAL_SRFLUX,'
#endif
#ifdef DJ_GRADPS
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
#endif
#ifdef DOUBLE_PRECISION
      IF (Master) WRITE (stdout,20) 'DOUBLE_PRECISION',                 &
     &   'Double precision arithmetic.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
#endif
#if defined DRENNAN
      IF (Master) WRITE (stdout,20) 'DRENNAN',                          &
     &   'Drennan et al (2003) relation for ZoW in bulk fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DRENNAN,'
#endif
#ifdef ECOSIM
      IF (Master) WRITE (stdout,20) 'ECOSIM',                           &
     &   'Bio-Optical EcoSim model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ECOSIM,'
      ibiology=ibiology+1
#endif
#if defined EMINUSP && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'EMINUSP',                          &
     &   'Compute Salt Flux using E-P.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' EMINUSP,'
#endif
#if defined ESMF_LIB && defined COAWST_COUPLING
      IF (Master) WRITE (stdout,20) 'ESMF_LIB',                         &
     &   'Using Earth System Modeling Framework library.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ESMF_LIB,'
#endif
#ifdef ENKF_RESTART
      IF (Master) WRITE (stdout,20) 'ENKF_RESTART',                     &
     &   'Writting fields for offline Ensemble Kalman Filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ENKF_RESTART,'
#endif
#ifdef ENSEMBLE
      IF (Master) WRITE (stdout,20) 'ENSEMBLE',                         &
     &   'Ensemble Forecasting Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ENSEMBLE,'
#endif
#if defined EVOLVED_LCZ && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'EVOLVED_LCZ',                      &
     &   'Compute 4D-Var evolved Hessian singular vectors.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' EVOLVED_LCZ,'
#endif
#ifdef FASTICE_CLIMATOLOGY
      IF (Master) WRITE(stdout,20) 'FASTICE_CLIMATOLOGY',               &
     &   'Climatology of landfast ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' FASTICE_CLIMATOLOGY,'
#endif
#ifdef FILTERED
      IF (Master) WRITE(stdout,20) 'FILTERED',                          &
     &   'Writing out time-filtered output.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' FILTERED,'
#endif
#ifdef FLOATS
      IF (Master) WRITE (stdout,20) 'FLOATS',                           &
     &   'Simulated Lagrangian drifters.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' FLOATS,'
#endif
#if defined FLOAT_BIOLOGY && defined FLOATS
      IF (Master) WRITE (stdout,20) 'FLOAT_BIOLOGY',                   &
     &   'Biological behavior on Lagrangian drifters.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FLOAT_BIOLOGY,'
#endif
#if defined FLOAT_OYSTER && defined FLOATS
      IF (Master) WRITE (stdout,20) 'FLOAT_OYSTER',                     &
     &   'Oyster model behavior on Lagrangian drifters.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' FLOAT_OYSTER,'
#endif
#if defined FLOAT_STICKY && defined FLOATS && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'FLOAT_STICKY',                     &
     &   'Reflect (stick) floats that hit the surface (bottom).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' FLOAT_STICKY,'
#endif
#if defined FLOATS && defined FLOAT_VWALK && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'FLOAT_VWALK',                      &
     &   'Lagrangian drifters with vertical diffusion, random walk.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' FLOAT_VWALK,'
#endif
#ifdef FORCING_SV
      IF (Master) WRITE (stdout,20) 'FORCING_SV',                       &
     &   'Forcing Singular Vectors Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' FORCING_SV,'
      idriver=idriver+1
#endif
#if defined FORWARD_MIXING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'FORWARD_MIXING',                   &
     &   'Read in Forward vertical mixing for Tangent/Adjoint.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' FORWARD_MIXING,'
#endif
#ifdef FORWARD_READ
      IF (Master) WRITE (stdout,20) 'FORWARD_READ',                     &
     &   'Read in Forward solution for Tangent/Adjoint.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' FORWARD_READ,'
#endif
#ifdef FORWARD_RHS
      IF (Master) WRITE (stdout,20) 'FORWARD_RHS',                      &
     &   'Process Forward RHS terms for Tangent/Adjoint.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' FORWARD_RHS,'
#endif
#ifdef FORWARD_WRITE
      IF (Master) WRITE (stdout,20) 'FORWARD_WRITE',                    &
     &   'Write out Forward solution for Tangent/Adjoint.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FORWARD_WRITE,'
#endif
#ifdef FSOBC_REDUCED
      IF (Master) WRITE (stdout,20) 'FSOBC_REDUCED',                    &
     &   'Using free-surface data in reduced physics conditions'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FSOBC_REDUCED,'
#endif
#ifdef FT_EIGENMODES
      IF (Master) WRITE (stdout,20) 'FT_EIGENMODES',                    &
     &   'Finite Time Eigenmodes: Normal Modes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FT_EIGENMODES,'
      idriver=idriver+1
#endif
#if defined FOUR_DVAR || defined PROPAGATOR || defined VERIFICATION
# ifdef FULL_GRID
      IF (Master) WRITE (stdout,20) 'FULL_GRID',                        &
#  if (defined FOUR_DVAR || defined VERIFICATION)
     &   'Considering observations at interior and boundary points.'
#  elif defined PROPAGATOR
     &   'State vector includes interior and boundary points.'
#  endif
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' FULL_GRID,'
# else
      IF (Master) WRITE (stdout,20) '!FULL_GRID',                       &
#  if (defined FOUR_DVAR || defined VERIFICATION)
     &   'Considering observations at interior points only.'
#  elif defined PROPAGATOR
     &   'State vector includes interior points only.'
#  endif
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' !FULL_GRID,'
# endif
#endif
#if defined GEOPOTENTIAL_HCONV && defined FOUR_DVAR && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'GEOPOTENTIAL_HCONV',               &
     &   'Horizontal convolutions along geopetential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' GEOPOTENTIAL_HCONV,'
#endif
#ifdef GLOBAL_PERIODIC
      IF (Master) WRITE (stdout,20) 'GLOBAL_PERIODIC',                  &
     &   'Globally periodic domain when interpolating forcing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' GLOBAL_PERIODIC,'
#endif
#ifdef GLS_MIXING
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef HDF5
      IF (Master) WRITE (stdout,20) 'HDF5',                             &
     &   'Creating NetCDF-4/HDF5 format files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' HDF5,'
#endif
#ifdef HESSIAN_FSV
      IF (Master) WRITE (stdout,20) 'HESSIAN_FSV',                      &
     &   'Hessian (4D-Var) Forcing Singular Vectors Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' HESSIAN_FSV,'
#endif
#ifdef HESSIAN_SO
      IF (Master) WRITE (stdout,20) 'HESSIAN_SO',                       &
     &   'Hessian (4D-Var) Stochastic Optimals Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' HESSIAN_SO,'
#endif
#ifdef HESSIAN_SV
      IF (Master) WRITE (stdout,20) 'HESSIAN_SV',                       &
     &   'Hessian (4D-Var) Singular Vectors Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' HESSIAN_SV,'
#endif
# ifdef HISTORY2
      IF (Master) WRITE (stdout,20) 'HISTORY2',                        &
     &   'Writing out secondary surface fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' HISTORY2,'
# endif
#if defined HOLLING_GRAZING && defined NEMURO
      IF (Master) WRITE (stdout,20) 'HOLLING_GRAZING',                  &
     &   'Holling-type s-shaped grazing algorithm, Implicit.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' HOLLING_GRAZING,'
#endif
#ifdef ICE_ADVECT
      IF (Master) WRITE (stdout,20) 'ICE_ADVECT',                       &
     &   'Advection of ice tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_ADVECT,'
#endif
#ifdef ICE_ALB_EC92
      IF (Master) WRITE (stdout,20) 'ICE_ALB_EC92',                     &
     &   'Ebert and Curry albedo formula.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_ALB_EC92,'
#endif
#ifdef ICE_BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'ICE_BULK_FLUXES',                  &
     &   'Ice bulk fluxes from the atmosphere.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ICE_BULK_FLUXES,'
#endif
#ifdef ICE_CONVSNOW
      IF (Master) WRITE (stdout,20) 'ICE_CONVSNOW',                     &
     &   'Conversion of flooded snow to ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_CONVSNOW,'
#endif
#ifdef ICE_EVP
      IF (Master) WRITE (stdout,20) 'ICE_EVP',                          &
     &   'Elastic-viscous-plastic ice rheology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ICE_EVP,'
#endif
#ifdef ICE_I_O
      IF (Master) WRITE (stdout,20) 'ICE_I_O',                          &
     &   'Put shortwave heat into ice interior.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ICE_I_O,'
#endif
#ifdef ICE_LANDFAST
      IF (Master) WRITE (stdout,20) 'ICE_LANDFAST',                    &
     &   'Lemieux shorefast ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_LANDFAST,'
#endif
#ifdef ICE_MK
      IF (Master) WRITE (stdout,20) 'ICE_MK',                           &
     &   'Mellor-Kantha ice thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ICE_MK,'
#endif
#ifdef ICE_MODEL
      IF (Master) WRITE (stdout,20) 'ICE_MODEL',                        &
     &   'Include sea ice model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ICE_MODEL,'
#endif
#ifdef ICE_MOMENTUM
      IF (Master) WRITE (stdout,20) 'ICE_MOMENTUM',                     &
     &   'Compute ice momentum equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_MOMENTUM,'
#endif
#ifdef ICE_SHOREFAST
      IF (Master) WRITE (stdout,20) 'ICE_SHOREFAST',                    &
     &   'Motionless shorefast ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ICE_SHOREFAST,'
#endif
#ifdef ICE_SMOLAR
      IF (Master) WRITE (stdout,20) 'ICE_SMOLAR',                       &
     &   'Advect ice tracers with MPDATA scheme.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_SMOLAR,'
#endif
#ifdef ICE_STRENGTH_QUAD
      IF (Master) WRITE (stdout,20) 'ICE_STRENGTH_QUAD',                &
     &   'Ice strength a quadratic function of thickness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' ICE_STRENGTH_QUAD,'
#endif
#ifdef ICE_THERMO
      IF (Master) WRITE (stdout,20) 'ICE_THERMO',                       &
     &   'Include ice thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_THERMO,'
#endif
#ifdef ICESHELF
      IF (Master) WRITE (stdout,20) 'ICESHELF',                         &
     &   'Include Ice Shelf Cavities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ICESHELF,'
#endif
#ifdef ICESHELF_3EQ
      IF (Master) WRITE (stdout,20) 'ICESHELF_3EQ',                     &
     &   'Include 3eq Ice Shelf Thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICESHELF_3EQ,'
#endif
#if defined IMPLICIT_VCONV && defined FOUR_DVAR && defined VCONVOLUTION
      IF (Master) WRITE (stdout,20) 'IMPLICIT_VCONV',                   &
     &   'Implicit Vertical Convolution Algorithm.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' IMPLICIT_VCONV,'
#endif
#ifdef IMPULSE
      IF (Master) WRITE (stdout,20) 'IMPULSE',                          &
     &   'Processing Adjoint Impulse forcing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' IMPULSE,'
#endif
#if defined LIMIT_STFLX_COOLING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'LIMIT_STFLX_COOLING',              &
     &   'Suppress further cooling if SST is at freezing point.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' LIMIT_STFLX_COOLING,'
#endif
#ifdef INNER_PRODUCT
      IF (Master) WRITE (stdout,20) 'INNER_PRODUCT',                    &
     &   'Tangent/Adjoint State Matrices Inner Product Test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' INNER_PRODUCT,'
      idriver=idriver+1
#endif
#if defined IRON_LIMIT && (NPZD_IRON || defined NEMURO || \
             defined BIO_GOANPZ || defined BIO_UMAINE)
      IF (Master) WRITE (stdout,20) 'IRON LIMITATION',                  &
     &   'Iron limitation on phytoplankton growth.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' IRON_LIMIT,'
#endif
#if defined IRON_RELAX && (defined NPZD_IRON || defined NEMURO || \
          defined BIO_UMAINE)
      IF (Master) WRITE (stdout,20) 'IRON_RELAX',                       &
     &   'Nudge dissolved iron to climatology'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' IRON_RELAX,'
#endif
#if defined IRON_RSIN && (defined NEMURO || defined BIO_UMAINE)
      IF (Master) WRITE (stdout,20) 'IRON_RSIN',                       &
     &   'Variable Si:N ratio based on iron'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' IRON_RSIN,'
#endif
#ifdef IS4DVAR
      IF (Master) WRITE (stdout,20) 'IS4DVAR',                          &
     &   'Incremental strong constraint 4D-Var data assimilation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' IS4DVAR,'
      idriver=idriver+1
#endif
#ifdef IS4DVAR_SENSITIVITY
      IF (Master) WRITE (stdout,20) 'IS4DVAR_SENSITIVITY',              &
     &   'I4D-Var Observations Sensitivity Analysis.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' IS4DVAR_SENSITIVITY,'
      idriver=idriver+1
#endif
#if defined INLINE_2DIO && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'INLINE_2DIO',                      &
     &   'Processing 3D IO level by level to reduce memory needs.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' INLINE_2DIO,'
#endif
#if defined IVLEV_EXPLICIT && defined NEMURO
      IF (Master) WRITE (stdout,20) 'IVLEV_EXPLICIT',                   &
     &   'Ivlev Grazing algorithm, Explicit.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' IVLEV_EXPLICIT,'
#endif
#if defined KANTHA_CLAYSON && (defined GLS_MIXING || defined MY25_MIXING)
      IF (Master) WRITE (stdout,20) 'KANTHA_CLAYSON',                   &
     &   'Kantha and Clayson stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' KANTHA_CLAYSON,'
#endif
#if defined LCZ_FINAL && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'LCZ_FINAL',                        &
     &   'Compute 4D-Var Hessian singular vectors.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LCZ_FINAL,'
#endif
#ifdef LMD_BKPP
      IF (Master) WRITE (stdout,20) 'LMD_BKPP',                         &
     &   'KPP bottom boundary layer mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_BKPP,'
#endif
#ifdef LMD_CONVEC
      IF (Master) WRITE (stdout,20) 'LMD_CONVEC',                       &
     &   'LMD convective mixing due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_CONVEC,'
#endif
#ifdef LMD_DDMIX
      IF (Master) WRITE (stdout,20) 'LMD_DDMIX',                        &
     &   'LMD double-diffusive mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LMD_DDMIX,'
#endif
#ifdef LMD_MIXING
      IF (Master) WRITE (stdout,20) 'LMD_MIXING',                       &
     &   'Large/McWilliams/Doney interior mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef LMD_NONLOCAL
      IF (Master) WRITE (stdout,20) 'LMD_NONLOCAL',                     &
     &   'LMD convective nonlocal transport.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LMD_NONLOCAL,'
#endif
#ifdef LMD_RIMIX
      IF (Master) WRITE (stdout,20) 'LMD_RIMIX',                        &
     &   'LMD diffusivity due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LMD_RIMIX,'
#endif
#ifdef LMD_SHAPIRO
      IF (Master) WRITE (stdout,20) 'LMD_SHAPIRO',                      &
     &   'Shapiro filtering boundary layer depth.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' LMD_SHAPIRO,'
#endif
#ifdef LMD_SKPP
      IF (Master) WRITE (stdout,20) 'LMD_SKPP',                         &
     &   'KPP surface boundary layer mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_SKPP,'
#endif
#if defined LONGWAVE && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'LONGWAVE',                         &
     &   'Compute net longwave radiation internally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LONGWAVE,'
#endif
#if defined LONGWAVE_OUT && defined BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'LONGWAVE_OUT',                     &
     &   'Compute outgoing longwave radiation internally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LONGWAVE_OUT,'
#endif
#if defined MIX_ISO_TS && (defined TS_DIF2 || defined TS_DIF4)
# ifdef MAX_SLOPE
      IF (Master) WRITE (stdout,20) 'MAX_SLOPE',                        &
     &   'Maximum density slope in isopycnic diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MAX_SLOPE,'
# endif
# ifdef MIN_STRAT
      IF (Master) WRITE (stdout,20) 'MIN_STRAT',                        &
     &   'Minimum density stratification in isopycnic diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MIN_STRAT,'
# endif
#endif
#ifdef M2TIDE_DIFF
      IF (Master) WRITE (stdout,20) 'M2TIDE_DIFF',                      &
     &   'Add barotropic tidal mixing parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' M2TIDE_DIFF,'
#endif
#ifdef MASKING
      IF (Master) WRITE (stdout,20) 'MASKING',                          &
     &   'Land/Sea masking.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
#endif
#ifdef MB_BBL
      IF (Master) WRITE (stdout,20) 'MB_BBL',                           &
     &   'Blaas Bottom Boundary Layer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' MB_BBL,'
      ibbl=ibbl+1
# ifdef MB_CALC_UB
      IF (Master) WRITE (stdout,20) 'MB_CALC_UB',                       &
     &   'Internal computation of bottom orbital velocity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MB_CALC_UB,'
# endif
# ifdef MB_CALC_ZNOT
      IF (Master) WRITE (stdout,20) 'MB_CALC_ZNOT',                     &
     &   'Internal computation of bottom roughness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' MB_CALC_ZNOT,'
# endif
# ifdef MB_Z0BIO
      IF (Master) WRITE (stdout,20) 'MB_Z0BIO',                         &
     &   'Biogenic bedform roughness ripple height and length.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MB_Z0BIO,'
# endif
# ifdef MB_Z0BL
      IF (Master) WRITE (stdout,20) 'MB_Z0BL',                          &
     &   'Bedload roughness for ripples.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MB_Z0BL,'
# endif
# ifdef MB_Z0RIP
      IF (Master) WRITE (stdout,20) 'MB_Z0RIP',                         &
     &   'Bedform roughness ripple height and length.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MB_Z0RIP,'
# endif
#endif
#if defined MCT_LIB && defined COAWST_COUPLING
      IF (Master) WRITE (stdout,20) 'MCT_LIB',                          &
     &   'Using Model Coupling Toolkit library.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MCT_LIB,'
#endif
#ifdef MELT_PONDS
      IF (Master) WRITE (stdout,20) 'MELT_PONDS',                       &
     &   'Melt ponds on ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MELT_PONDS,'
#endif
#ifdef MICLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'MICLIMATOLOGY',                    &
     &   'Ice momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' MICLIMATOLOGY,'
#endif
#ifdef MICLM_NUDGING
      IF (Master) WRITE (stdout,20) 'MICLM_NUDGING',                    &
     &   'Nudging to ice momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' MICLM_NUDGING,'
#endif
#if defined MINRES && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'MINRES',                           &
     &   'Minimal Residual Method for 4D-Var minimization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' MINRES,'
#endif
#if (defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D
# ifdef MIX_GEO_TS
      IF (Master) WRITE (stdout,20) 'MIX_GEO_TS',                       &
     &   'Mixing of tracers along geopotential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_GEO_TS,'
# endif
# ifdef MIX_ISO_TS
      IF (Master) WRITE (stdout,20) 'MIX_ISO_TS',                       &
     &   'Mixing of tracers along isopycnal surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_ISO_TS,'
# endif
# ifdef MIX_S_TS
      IF (Master) WRITE (stdout,20) 'MIX_S_TS',                         &
     &   'Mixing of tracers along constant S-surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_TS,'
# endif
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && defined SOLVE3D
# ifdef MIX_GEO_UV
      IF (Master) WRITE (stdout,20) 'MIX_GEO_UV',                       &
     &   'Mixing of momentum along geopotential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MIX_GEO_UV,'
# endif
# ifdef MIX_S_UV
      IF (Master) WRITE (stdout,20) 'MIX_S_UV',                         &
     &   'Mixing of momentum along constant S-surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_UV,'
# endif
#endif
#ifdef MPI
      IF (Master) WRITE (stdout,20) 'MPI',                              &
     &   'MPI distributed-memory configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
#endif
#if defined MULTIPLE_TLM && defined TANGENT
      IF (Master) WRITE (stdout,20) 'MULTIPLE_TLM',                     &
     &   'Create multiple TLM history files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' MULTIPLE_TLM,'
#endif
#ifdef MY25_MIXING
      IF (Master) WRITE (stdout,20) 'MY25_MIXING',                      &
     &   'Mellor/Yamada Level-2.5 mixing closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' MY25_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef WEC_MELLOR
      IF (Master) WRITE (stdout,20) 'WEC_MELLOR',                       &
     &   'Mellore wave current interaction terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' WEC_MELLOR,'
      nearshore=nearshore+1
#endif
#ifdef WEC_VF
      IF (Master) WRITE (stdout,20) 'WEC_VF',                           &
     &   'Vortex Force wave current interaction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' WEC_VF,'
      nearshore=nearshore+1
#endif
#ifdef BOTTOM_STREAMING
      IF (Master) WRITE (stdout,20) 'BOTTOM_STREAMING',                 &
     &   'Wave current streaming term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' BOTTOM_STREAMING,'
#endif
#ifdef SURFACE_STREAMING
      IF (Master) WRITE (stdout,20) 'SURFACE_STREAMING',                 &
     &   'Wave current streaming term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SURFACE_STREAMING,'
#endif
#ifdef WAVE_MIXING
      IF (Master) WRITE (stdout,20) 'WAVE_MIXING',                      &
     &   'Enhanced veritcal mixing from waves term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' WEC_WAVE_MIX,'
#endif
#ifdef WDISS_THORGUZA
      IF (Master) WRITE (stdout,20) 'WDISS_THORGUZA',                   &
     &   'Wave energy dissipation computed from Thorton/Guza 1986.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WDISS_THORGUZA,'
#endif
#ifdef WDISS_CHURTHOR
      IF (Master) WRITE (stdout,20) 'WDISS_CHURTHOR',                   &
     &   'Wave energy dissipation computed from Church/Thorton 1986.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WDISS_CHURTHOR,'
#endif
#ifdef WDISS_WAVEMOD
      IF (Master) WRITE (stdout,20) 'WDISS_WAVEMOD',                    &
     &   'Wave energy dissipation acquired from coupled wave model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' WDISS_WAVEMOD,'
#endif
#ifdef ROLLER_SVENDSEN
      IF (Master) WRITE (stdout,20) 'ROLLER_SVENDSEN',                  &
     &   'Wave energy roller based on Svendsen 1984.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ROLLER_SVENDSEN,'
#endif
#ifdef ROLLER_MONO
      IF (Master) WRITE (stdout,20) 'ROLLER_MONO',                      &
     &   'Wave energy roller from monochromatic waves.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ROLLER_MONO,'
#endif
#ifdef ROLLER_RENIERS
      IF (Master) WRITE (stdout,20) 'ROLLER_RENIERS',                   &
     &   'Wave energy roller based on Reniers 2004.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ROLLER_RENIERS,'
#endif
#ifdef NEMURO
      IF (Master) WRITE (stdout,20) 'NEMURO',                           &
     &   'Nemuro Lower Trophic Level Ecosystem Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' NEMURO,'
      ibiology=ibiology+1
#endif
#ifdef NEMURO_SED1
      IF (Master) WRITE (stdout,20) 'NEMURO_SED1',                      &
     &   'Nemuro with sediment remineralization type 1.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NEMURO_SED1,'
#endif
#ifdef NESTING
      IF (Master) WRITE (stdout,20) 'NESTING',                          &
     &   'Nesting grids: Composite and Refinement.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' NESTING,'
#endif
#if defined NL_BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'NL_BULK_FLUXES',                   &
     &   'Using bulk fluxes computed by nonlinear model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' NL_BULK_FLUXES,'
#endif
#if defined NLM_OUTER && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'NLM_OUTER',                        &
     &   'Using the Nonlinear model as Basic State in Outer Loop.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' NLM_OUTER,'
#endif
#ifdef NONLINEAR
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
#endif
#ifdef SOLVE3D
# if defined NONLIN_EOS
      IF (Master) WRITE (stdout,20) 'NONLIN_EOS',                       &
     &   'Nonlinear Equation of State for seawater.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NONLIN_EOS,'
# else
      IF (Master) WRITE (stdout,20) '!NONLIN_EOS',                      &
     &   'Linear Equation of State for seawater.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !NONLIN_EOS,'
# endif
#endif
#if defined NO_CORRECT_TRACER && defined NESTING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'NO_CORRECT_TRACER',                &
     &   'Avoid two-way correction of coaser grid boundary tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NO_CORRECT_TRACER,'
#endif
#ifdef NO_HIS
      IF (Master) WRITE (stdout,20) 'NO_HIS',                           &
     &   'Turn off writing of history files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' NO_HIS,'
#endif
#ifdef NO_LBC_ATT
      IF (Master) WRITE (stdout,20) 'NO_LBC_ATT',                       &
     &   'Not checking NetCDF global attribute NLM_LBC during restart.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NO_LBC_ATT,'
#endif
#if defined NO_READ_GHOST && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'NO_READ_GHOST',                    &
     &   'Not including ghost points during scattering read data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NO_READ_GHOST,'
#endif
#if defined NO_SCORRECTION_ICE
      IF (Master) WRITE (stdout,20) 'NO_SCORRECTION_ICE',               &
     &   'Not including SCORRECTION under the ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' NO_SCORRECTION_ICE,'
#endif
#ifdef NO_SNOW
      IF (Master) WRITE (stdout,20) 'NO_SNOW',                           &
     &   'Turn off snow in sea-ice model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' NO_SNOW,'
#endif
#ifdef NO_WRITE_GRID
      IF (Master) WRITE (stdout,20) 'NO_WRITE_GRID',                    &
     &   'Not Writing grid arrays into NetCDF ouput files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NO_WRITE_GRID,'
#endif
#ifdef NPZD_FRANKS
      IF (Master) WRITE (stdout,20) 'NPZD_FRANKS',                      &
     &   'NPZD Biological Model, Franks et al. fomulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NPZD_FRANKS,'
      ibiology=ibiology+1
#endif
#ifdef NPZD_IRON
      IF (Master) WRITE (stdout,20) 'NPZD_IRON',                        &
     &   'NPZD Biological Model (Powell et al.) with iron limitation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' NPZD_IRON,'
      ibiology=ibiology+1
#endif
#ifdef NPZD_POWELL
      IF (Master) WRITE (stdout,20) 'NPZD_POWELL',                      &
     &   'NPZD Biological Model, Powell et al. fomulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NPZD_POWELL,'
      ibiology=ibiology+1
#endif
#if defined N2S2_HORAVG && (defined GLS_MIXING || defined MY25_MIXING)
      IF (Master) WRITE (stdout,20) 'N2S2_HORAVG',                      &
     &   'Horizontal smoothing of buoyancy and shear.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
#endif
#ifdef OBSERVATIONS
      IF (Master) WRITE (stdout,20) 'OBSERVATIONS',                     &
     &   'Processing 4DVar observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' OBSERVATIONS,'
#endif
#ifdef OBS_IMPACT
      IF (Master) WRITE (stdout,20) 'OBS_IMPACT',                       &
     &   'Compute observation impact to the 4DVAR system.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' OBS_IMPACT,'
#endif
#if defined OBS_IMPACT && defined OBS_IMPACT_SPLIT
      IF (Master) WRITE (stdout,20) 'OBS_IMPACT_SPLIT',                 &
     &   'Separate impact due to IC, surface forcing, and OBC'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' OBS_IMPACT_SPLIT,'
#endif
#ifdef OCLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'OCLIMATOLOGY',                     &
     &   'Processing omega climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' OCLIMATOLOGY,'
#endif
#if defined OFFLINE_FLOATS && defined CLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'OFFLINE_FLOATS',                   &
     &   'Offline run (use model fields from history files).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' OFFLINE_FLOATS,'
#endif
#if defined OFFLINE_BIOLOGY && defined CLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'OFFLINE_BIOLOGY',                  &
     &   'Offline run (use model fields from history files).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' OFFLINE_BIOLOGY,'
#endif
#if defined ONE_TRACER_SOURCE
      IF (Master) WRITE (stdout,20) 'ONE_TRACER_SOURCE',                &
     &   'One value per tracer for all sources.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' ONE_TRACER_SOURCE,'
#endif
#ifdef NESTING
# ifdef ONE_WAY
      IF (Master) WRITE (stdout,20) 'ONE_WAY',                          &
     &   'One-way nesting in refinement grids.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ONE_WAY,'
# else
      IF (Master) WRITE (stdout,20) '!ONE_WAY',                         &
     &   'Two-way nesting in refinement grids.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' !ONE_WAY,'
# endif
#endif
#ifdef OUT_DOUBLE
      IF (Master) WRITE (stdout,20) 'OUT_DOUBLE',                       &
     &   'Double precision output fields in NetCDF files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' OUT_DOUBLE,'
#endif
#ifdef OUTFLOW_MASK
      IF (Master) WRITE (stdout,20) 'OUTFLOW_MASK',                     &
     &   'Hibler-style outflow ice cells.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' OUTFLOW_MASK,'
#endif
#ifdef _OPENMP
      IF (Master) WRITE (stdout,20) '_OPENMP',                          &
     &   'OpenMP parallel shared-memory directives.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' _OPENMP,'
#endif
#ifdef OPT_OBSERVATIONS
      IF (Master) WRITE (stdout,20) 'OPT_OBSERVATIONS',                 &
     &   'Optimal Observations Driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' OPT_OBSERVATIONS,'
      idriver=idriver+1
#endif
#ifdef OPT_PERTURBATION
      IF (Master) WRITE (stdout,20) 'OPT_PERTURBATION',                 &
     &   'Optimal Perturbation Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' OPT_PERTURBATION,'
      idriver=idriver+1
#endif
#if defined OXYGEN && defined BIO_FENNEL
# ifdef OCMIP_OXYGEN_SC
      IF (Master) WRITE (stdout,20) 'OCMIP_OXYGEN_SC',                  &
     &   'Schmidt number from Keeling et al. (1998).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' OCMIP_OXYGEN_SC,'
# else
      IF (Master) WRITE (stdout,20) '!OCMIP_OXYGEN_SC',                 &
     &   'Schmidt number from Wanninkhof (1992).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !OCMIP_OXYGEN_SC,'
# endif
      IF (Master) WRITE (stdout,20) 'OXYGEN',                           &
     &   'Add oxygen dynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' OXYGEN,'
#endif
#if defined OXYGEN && defined ESTUARYBGC 
# ifdef OCMIP_OXYGEN_SC
      IF (Master) WRITE (stdout,20) 'OCMIP_OXYGEN_SC',                  &
     &   'Schmidt number from Keeling et al. (1998).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' OCMIP_OXYGEN_SC,'
# else
      IF (Master) WRITE (stdout,20) '!OCMIP_OXYGEN_SC',                 &
     &   'Schmidt number from Wanninkhof (1992).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !OCMIP_OXYGEN_SC,'
# endif
      IF (Master) WRITE (stdout,20) 'OXYGEN',                           &
     &   'Add oxygen dynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' OXYGEN,'
#endif
#if defined PARALLEL_IN && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'PARALLEL_IN',                      &
     &   'Parallel I/O processing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' PARALLEL_IN,'
#endif
#if defined PARALLEL_OUT && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'PARALLEL_OUT',                      &
     &   'Parallel I/O processing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' PARALLEL_OUT,'
#endif
#if defined PARALLEL_IO && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'PARALLEL_IO',                      &
     &   'Parallel I/O processing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' PARALLEL_IO,'
#endif
#ifdef PERFECT_RESTART
      IF (Master) WRITE (stdout,20) 'PERFECT_RESTART',                  &
     &   'Processing perfect restart variables.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' PERFECT_RESTART,'
#endif
#ifdef PICARD_TEST
      IF (Master) WRITE (stdout,20) 'PICARD_TEST',                      &
     &   'Representers tangent linear model Picard iterations test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' PICARD TEST,'
      idriver=idriver+1
#endif
#ifdef PJ_GRADP
      IF (Master) WRITE (stdout,20) 'PJ_GRADP',                         &
     &   'Finite volume Pressure Jacobian (Lin, 1997).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' PJ_GRADP,'
#endif
#ifdef PJ_GRADPQ2
      IF (Master) WRITE (stdout,20) 'PJ_GRADPQ2',                       &
     &   'Quartic 2 polynomial Pressure Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' PJ_GRADPQ2,'
#endif
#ifdef PJ_GRADPQ4
      IF (Master) WRITE (stdout,20) 'PJ_GRADPQ4',                       &
     &   'Quartic 4 polynomial Pressure Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' PJ_GRADPQ4,'
#endif
#ifdef POSITIVE_ZERO
      IF (Master) WRITE (stdout,20) 'POSITIVE_ZERO',                    &
     &   'Impose positive zero in output data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' POSITIVE_ZERO,'
#endif
#if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'POSTERIOR_EOFS',                   &
     &   'Estimating posterior analysis error covariace matrix EOFs.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' POSTERIOR_EOFS,'
#endif
#if defined POSTERIOR_ERROR_F && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'POSTERIOR_ERROR_F',                &
     &   'Estimating final posterior analysis error covariace matrix.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' POSTERIOR_ERROR_F,'
#endif
#if defined POSTERIOR_ERROR_I && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'POSTERIOR_ERROR_I',                &
     &   'Estimating initial posterior analysis error covariace matrix.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' POSTERIOR_ERROR_F,'
#endif
#ifdef POT_TIDES
      IF (Master) WRITE (stdout,20) 'POT_TIDES',                        &
     &   'Add tidal potential terms to the forcing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POT_TIDES,'
#endif
#if defined POWER_LAW && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
#endif
#if !(defined PJ_GRADPQ4 || defined PJ_GRADPQ2 || defined PJ_GRADP || \
      defined DJ_GRADPS) && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'PRSGRD31',                         &
     &   'Standard density Jacobian formulation (Song, 1998).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' PRSGRD31,'
#endif
#ifdef PRIMARY_PROD
      IF (Master) WRITE (stdout,20) 'PRIMARY_PROD',                      &
     &   'Nemuro or CoSiNE Net primary production.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' PRIMARY_PROD,'
#endif
#ifdef PROFILE
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
#endif
#ifdef PSEUDOSPECTRA
      IF (Master) WRITE (stdout,20) 'PSEUDOSPECTRA',                    &
     &   'Pseudospectra Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' PSEUDOSPECTRA,'
#endif
#ifdef QCORRECTION
      IF (Master) WRITE (stdout,20) 'QCORRECTION',                      &
     &   'Surface net heat flux correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' QCORRECTION,'
#endif
#if defined RUNOFF
      IF (Master) WRITE (stdout,20) 'RUNOFF',                           &
     &   'Surface flux of ungauged river runoff.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' RUNOFF,'
#endif
#if defined RUNOFF_SSH
      IF (Master) WRITE (stdout,20) 'RUNOFF_SSH',                       &
     &   'Surface flux of ungauged river runoff affecting SSH.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RUNOFF_SSH,'
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# if defined K_C2ADVECTION
      IF (Master) WRITE (stdout,20) 'K_C2ADVECTION',                    &
     &   'Second-order centered differences advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' K_C2ADVECTION,'
# elif defined K_C4ADVECTION
      IF (Master) WRITE (stdout,20) 'K_C4ADVECTION',                    &
     &   'Fourth-order centered differences advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' K_C4ADVECTION,'
# else
      IF (Master) WRITE (stdout,20) 'K_GSCHEME',                        &
     &   'Third-order upstream advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' K_GSCHEME,'
# endif
# ifdef TKE_DIF2
      IF (Master) WRITE (stdout,20) 'TKE_DIF2',                         &
     &   'Harmonic mixing of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TKE_DIF2,'
# endif
# ifdef TKE_DIF4
      IF (Master) WRITE (stdout,20) 'TKE_DIF4',                         &
     &   'Biharmonic mixing of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TKE_DIF4,'
# endif
#endif
#ifdef RADIATION_2D
      IF (Master) WRITE (stdout,20) 'RADIATION_2D',                     &
     &   'Use tangential phase speed in radiation conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' RADIATION_2D,'
#endif
#if defined RAMP_TIDES && (defined SSH_TIDES || defined UV_TIDES)
      IF (Master) WRITE (stdout,20) 'RAMP_TIDES',                       &
     &   'Ramping tidal forcing for one day.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' RAMP_TIDES,'
#endif
#if defined READ_WATER && defined MASKING
      IF (Master) WRITE (stdout,20) 'READ_WATER',                       &
     &   'Reading data at water points only.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' READ_WATER,'
#endif
#if defined RECOMPUTE_4DVAR && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'RECOMPUTE_4DVAR',                  &
     &   'Recomputing 4DVar during data assimilation analysis.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' RECOMPUTE_4DVAR,'
#endif
#ifdef RED_TIDE
      IF (Master) WRITE (stdout,20) 'RED_TIDE',                         &
     &   'Red tide biological model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' RED_TIDE,'
      ibiology=ibiology+1
#endif
#ifdef R_SYMMETRY
      IF (Master) WRITE (stdout,20) 'REP_MATRIX',                       &
     &   'Representers Matrix Symmetry Test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' REP_MATRIX,'
      idriver=idriver+1
#endif
#if defined LMD_MIXING
# ifdef RI_HORAVG
      IF (Master) WRITE (stdout,20) 'RI_HORAVG',                        &
     &   'Smooth Richardson number horizontally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' RI_HORAVG,'
# endif
# ifdef RI_VERAVG
      IF (Master) WRITE (stdout,20) 'RI_VERAVG',                        &
     &   'Smooth Richardson number vertically.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' RI_VERAVG,'
# endif
#endif
#if defined SOLVE3D    && \
   (defined GLS_MIXING || defined LMD_MIXING || defined MY25_MIXING)
# ifdef RI_SPLINES
      IF (Master) WRITE (stdout,20) 'RI_SPLINES',                       &
     &   'Parabolic Spline Reconstruction for Richardson Number.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RI_SPLINES,'
# endif
#endif
#if !(defined PJ_GRADPQ4 || defined PJ_GRADPQ2 || defined PJ_GRADP || \
      defined DJ_GRADPS) && defined RHO_SURF   && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'RHO_SURF',                         &
     &   'Include difference between rho0 and surface density.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' RHO_SURF,'
#endif
#if defined RPCG && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'RPCG',                             &
     &   'Restricted B-preconditioned Lanczos minimization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' RPCG,'
#endif
#if defined RP_AVERAGES && defined TL_IOMS
      IF (Master) WRITE (stdout,20) 'RP_AVERAGES',                      &
     &   'Writing out time-averaged representer model fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' RP_AVERAGES,'
#endif
#ifdef RPM_DRIVER
      IF (Master) WRITE (stdout,20) 'RPM_DRIVER',                       &
     &   'Generic representers tangent linear model driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RPM_DRIVER,'
      idriver=idriver+1
#endif
#ifdef RPM_RELAXATION
      IF (Master) WRITE (stdout,20) 'RPM_RELAXATION',                   &
     &   'Diffusive Relaxation of RPM in Picard Iterations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' RPM_RELAXATION,'
#endif
#ifdef RST_SINGLE
      IF (Master) WRITE (stdout,20) 'RST_SINGLE',                       &
     &   'Single precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RST_SINGLE,'
#else
      IF (Master) WRITE (stdout,20) '!RST_SINGLE',                      &
     &   'Double precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
#endif
#ifdef SALINITY
      IF (Master) WRITE (stdout,20) 'SALINITY',                         &
     &   'Using salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
#endif
#ifdef SANITY_CHECK
      IF (Master) WRITE (stdout,20) 'SANITY_CHECK',                     &
     &   'Tangent Linear and Adjoint Models Correctness Check.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SANITY_CHECK,'
      idriver=idriver+1
#endif
#ifdef SCORRECTION
      IF (Master) WRITE (stdout,20) 'SCORRECTION',                      &
     &   'Surface salinity flux correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SCORRECTION,'
#endif
#ifdef SEDIMENT
      IF (Master) WRITE (stdout,20) 'SEDIMENT',                         &
     &   'Cohesive and noncohesive sediments.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SEDIMENT,'
# ifdef SED_DENS
      IF (Master) WRITE (stdout,20) 'SED_DENS',                         &
     &   'Include density of suspended sediment in equation of state.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SED_DENS,'
# endif
# ifdef SED_MORPH
      IF (Master) WRITE (stdout,20) 'SED_MORPH',                        &
     &   'Allow bottom model elevation to evolve.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SED_MORPH,'
# endif
# ifdef SUSPLOAD
      IF (Master) WRITE (stdout,20) 'SUSPLOAD',                         &
     &   'Activate suspended sediment transport.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SUSPLOAD,'
# endif
# ifdef SED_BIODIFF
      IF (Master) WRITE (stdout,20) 'SED_BIODIFF',                      &
     &   'Activate sediment bed biodiffusivity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SED_BIODIFF,'
# endif
# ifdef MIXED_BED
      IF (Master) WRITE (stdout,20) 'MIXED_BED',                        &
     &   'Activate mixed sediment bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MIXED_BED,'
      sbed_type=sbed_type+1
# endif
# ifdef COHESIVE_BED
      IF (Master) WRITE (stdout,20) 'COHESIVE_BED',                     &
     &   'Activate cohesive sediment bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' COHESIVE_BED,'
      sbed_type=sbed_type+1
# endif
# if defined NONCOHESIVE_BED1  || (!defined MIXED_BED && \
    !defined NONCOHESIVE_BED2  && !defined COHESIVE_BED)
      IF (Master) WRITE (stdout,20) 'NONCOHESIVE_BED1',                 &
     &   'Activate non-cohesive (original) sediment bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NONCOHESIVE_BED1,'
      sbed_type=sbed_type+1
# endif
# ifdef NONCOHESIVE_BED2
      IF (Master) WRITE (stdout,20) 'NONCOHESIVE_BED2',                 &
     &   'Activate non-cohesive (modified) sediment bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NONCOHESIVE_BED2,'
      sbed_type=sbed_type+1
# endif
# ifdef SED_FLOCS
      IF (Master) WRITE (stdout,20) 'SED_FLOCS',                        &
     &   'Activate sediment flocculation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SED_FLOCS,'
# endif
# ifdef FLOC_TURB_DISS
      IF (Master) WRITE (stdout,20) 'FLOC_TURB_DISS',                   &
     &   'Activate dissipation for flocculation based on turbulence.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' FLOC_TURB_DISS,'
# endif
# ifdef FLOC_BBL_DISS
      IF (Master) WRITE (stdout,20) 'FLOC_BBL_DISS',                    &
     &   'Activate dissipation for flocculation based on BBL.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FLOC_BBL_DISS,'
# endif
# ifdef SED_DEFLOC
      IF (Master) WRITE (stdout,20) 'SED_DEFLOC',                       &
     &   'Activate sediment deflocculation in bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' SED_DEFLOC,'
# endif
# ifdef SED_TAU_CD_CONST
      IF (Master) WRITE (stdout,20) 'SED_TAU_CD_CONST',                 &
     &   'Activate constant critical stress for deposition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SED_TAU_CD_CONST,'
# endif
# ifdef SED_TAU_CD_LIN
      IF (Master) WRITE (stdout,20) 'SED_TAU_CD_LIN',                   &
     &   'Activate linear critical stress for deposition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SED_TAU_CD_LIN,'
# endif
#endif
#ifdef SG_BBL
      IF (Master) WRITE (stdout,20) 'SG_BBL',                           &
     &   'Styles and Glenn Bottom Boundary Layer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' SG_BBL,'
      ibbl=ibbl+1
# ifdef SG_CALC_UB
      IF (Master) WRITE (stdout,20) 'SG_CALC_UB',                       &
     &   'Internal computation of bottom orbital velocity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' SG_CALC_UB,'
# endif
# ifdef SG_CALC_ZNOT
      IF (Master) WRITE (stdout,20) 'SG_CALC_ZNOT',                     &
     &   'Internal computation of bottom roughness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SG_CALC_ZNOT,'
# endif
# ifdef SG_LOGINT
      IF (Master) WRITE (stdout,20) 'SG_LOGINT',                        &
     &   'Bottom currents logarithmic interpolation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SG_LOGINT,'
# endif
#endif
#ifdef SNOWFALL
      IF (Master) WRITE (stdout,20) 'SNOWFALL',                         &
     &   'Surface snow accumulation on ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SNOWFALL,'
#endif
#ifdef SNOW_FROM_RAIN
      IF (Master) WRITE (stdout,20) 'SNOW_FROM_RAIN',                   &
     &   'Subtract snow from total precip (rain).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SNOW_FROM_RAIN,'
#endif
#ifdef SRELAXATION
      IF (Master) WRITE (stdout,20) 'SRELAXATION',                      &
     &   'Surface salinity relaxation as surface flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SRELAXATION,'
#endif
#ifdef SOLAR_SOURCE
      IF (Master) WRITE (stdout,20) 'SOLAR_SOURCE',                     &
     &   'Solar Radiation Source Term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SOLAR_SOURCE,'
#endif
#ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
#endif
#ifdef SO_SEMI
      IF (Master) WRITE (stdout,20) 'SO_SEMI',                          &
     &   'Stochastic Optimals Propagator, Semi-norm Estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SO_SEMI,'
      idriver=idriver+1
#endif
#ifdef SO_TRACE
      IF (Master) WRITE (stdout,20) 'SO_TRACE',                         &
     &   'Stochastic Optimals Propagator, Randomized Trace Estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SO_TRACE,'
      idriver=idriver+1
#endif
#ifdef SO_SEMI
# ifdef SO_SEMI_WHITE
      IF (Master) WRITE (stdout,20) 'SO_SEMI_WHITE',                    &
     &   'Stochastic Optimals, white noise processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SO_SEMI_WHITE,'
# else
      IF (Master) WRITE (stdout,20) '!SO_WHITE',                        &
     &   'Stochastic Optimals, red noise processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' !SO_SEMI_WHITE,'
# endif
#endif
#if defined SPLINES_VCONV && defined FOUR_DVAR && defined VCONVOLUTION
      IF (Master) WRITE (stdout,20) 'IMPLICIT_VCONV',                   &
     &   'Implicit, Parabolic Splines Vertical Convolution Algorithm.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SPLINES_VCONV,'
#endif
#if defined SPLINES_VDIFF && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'SPLINES_VDIFF',                    &
     &   'Parabolic Spline Reconstruction for Vertical Diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SPLINES_VDIFF,'
#endif
#if defined SPLINES_VVISC && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'SPLINES_VVISC',                    &
     &   'Parabolic Spline Reconstruction for Vertical Viscosity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SPLINES_VVISC,'
#endif
#ifdef SPHERICAL
      IF (Master) WRITE (stdout,20) 'SPHERICAL',                        &
     &   'Spherical grid configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SPHERICAL,'
#endif
#ifdef SRELAXATION
      IF (Master) WRITE (stdout,20) 'SRELAXATION',                      &
     &   'Surface salinity relaxation as surface flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SRELAXATION,'
#endif
#if defined SSH_TIDES
      IF (Master) WRITE (stdout,20) 'SSH_TIDES',                        &
     &   'Add tidal elevation to SSH climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SSH_TIDES,'
#endif
#if defined SSSC_THRESHOLD
      IF (Master) WRITE (stdout,20) 'SSSC_THRESHOLD',                   &
     &   'Limit to SCORRECTION nudging.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SSSC_THRESHOLD,'
#endif
#ifdef SSW_BBL
       IF (Master) WRITE (stdout,20) 'SSW_BBL',                         &
     &   'Styles and Glenn Bottom Boundary Layer - modified.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+8)=' SSW_BBL,'
       ibbl=ibbl+1
# ifdef SSW_CALC_UB
       IF (Master) WRITE (stdout,20) 'SSW_CALC_UB',                     &
     &   'Internal computation of bottom orbital velocity.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+14)=' SSW_CALC_UB,'
# endif
# ifdef SSW_CALC_ZNOT
       IF (Master) WRITE (stdout,20) 'SSW_CALC_ZNOT',                   &
     &   'Internal computation of bottom roughness.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+14)=' SSW_CALC_ZNOT,'
# endif
# ifdef SSW_FORM_DRAG_COR
       IF (Master) WRITE (stdout,20) 'SSW_FORM_DRAG_COR',               &
     &   'Form drag coefficient.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+19)=' SSW_FORM_DRAG_COR,'
# endif
# ifdef SSW_ZOBIO
       IF (Master) WRITE (stdout,20) 'SSW_Z0BIO',                       &
     &   'Biogenic bedfrom roughness from ripples.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+11)=' SSW_Z0BIO,'
# endif
# ifdef SSW_ZOBL
       IF (Master) WRITE (stdout,20) 'SSW_Z0BL',                        &
     &   'Bedload roughness from ripples.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+10)=' SSW_Z0BL,'
# endif
# ifdef SSW_ZORIP
       IF (Master) WRITE (stdout,20) 'SSW_Z0RIP',                       &
     &   'Bedfrom roughness from ripples.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+11)=' SSW_Z0RIP,'
# endif
# ifdef SSW_LOGINT
       IF (Master) WRITE (stdout,20) 'SSW_LOGINT',                      &
     &   'Bottom currents logarithmic interpolation.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+11)=' SSW_LOGINT,'
# endif
#endif
#ifdef STATIONS
      IF (Master) WRITE (stdout,20) 'STATIONS',                         &
     &   'Writing out station data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' STATIONS,'
#endif
#if defined STATIONS_CGRID && defined STATIONS
      IF (Master) WRITE (stdout,20) 'STATIONS_CGRID',                   &
     &   'Extracting station data at native C-grid locations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' STATIONS_CGRID,'
#endif
#ifdef STOCHASTIC_OPT
      IF (Master) WRITE (stdout,20) 'STOCHASTIC_OPT',                   &
     &   'Stochastic Optimals Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' STOCHASTIC_OPT,'
      idriver=idriver+1
# ifdef STOCH_OPT_WHITE
      IF (Master) WRITE (stdout,20) 'STOCH_OPT_WHITE',                  &
     &   'Stochastic Optimals, white noise processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' STOCH_OPT_WHITE,'
# else
      IF (Master) WRITE (stdout,20) '!STOCH_OPT_WHITE',                 &
     &   'Stochastic Optimals, red noise processes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !STOCH_OPT_WHITE,'
# endif
#endif
#ifdef SWAN_COUPLING
      IF (Master) WRITE (stdout,20) 'SWAN_COUPLING',                    &
     &   'SWAN model coupling.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SWAN_COUPLING,'
#endif
#if defined CARBON && defined BIO_FENNEL
# ifdef TALK_NONCONSERV
      IF (Master) WRITE (stdout,20) 'TALK_NONCONSERV',                  &
     &   'Alkalinity is affected by changes in nitrate or ammonium.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TALK_NONCONSERV,'
# else
      IF (Master) WRITE (stdout,20) '!TALK_NONCONSERV',                 &
     &   'Alkalinity is passive and unaffected by nitrate or ammonium.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !TALK_NONCONSERV,'
# endif
#endif
#if defined CARBON && defined ESTUARYBGC
# ifdef TALK_NONCONSERV
      IF (Master) WRITE (stdout,20) 'TALK_NONCONSERV',                  &
     &   'Alkalinity is affected by changes in nitrate or ammonium.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TALK_NONCONSERV,'
# else
      IF (Master) WRITE (stdout,20) '!TALK_NONCONSERV',                 &
     &   'Alkalinity is passive and unaffected by nitrate or ammonium.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !TALK_NONCONSERV,'
# endif
#endif
#ifdef TANGENT
      IF (Master) WRITE (stdout,20) 'TANGENT',                          &
     &   'Tangent Linear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TANGENT,'
#endif
#if defined TIME_CONV && defined WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'TIME_CONV',                        &
     &   'Weak-constraint 4D-Var time convolution.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' TIME_CONV,'
#endif
#if defined TIDES_ASTRO
      IF (Master) WRITE (stdout,20) 'TIDES_ASTRO',                      &
     &   'Compute astronomical (long-period tide) corrections.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' TIDES_ASTRO,'
#endif
#if defined TIME_INTERP_FLUX && defined NESTING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TIME_INTERP_FLUX',                 &
     &   'Time interpolate coaser grid mass flux in refinement .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' TIME_INTERP_FLUX,'
#endif
#if defined TL_AVERAGES && defined TANGENT
      IF (Master) WRITE (stdout,20) 'TL_AVERAGES',                      &
     &   'Writing out time-averaged tangent linear model fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' TL_AVERAGES,'
#endif
#ifdef TLM_DRIVER
      IF (Master) WRITE (stdout,20) 'TLM_DRIVER',                       &
     &   'Generic tangent linear model driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' TLM_DRIVER,'
      idriver=idriver+1
#endif
#if defined GLS_MIXING && defined TKE_WAVEDISS
      IF (Master) WRITE (stdout,20) 'TKE_WAVEDISS',                     &
     &   'Wave breaking surface tke flux based on wave amplitude.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TKE_WAVEDISS,'
#endif
#ifdef TL_IOMS
      IF (Master) WRITE (stdout,20) 'TL_IOMS',                          &
     &   'Representers Tangent Linear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TL_IOMS,'
#endif
#ifdef TLM_CHECK
      IF (Master) WRITE (stdout,20) 'TLM_CHECK',                        &
     &   'Tangent Linear Model Linearization Check.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' TLM_CHECK,'
#endif
#ifdef SOLVE3D
# if defined T_PASSIVE
      IF (Master) WRITE (stdout,20) 'T_PASSIVE',                        &
     &   'Advecting and diffusing inert passive tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' T_PASSIVE,'
# endif
# ifdef TS_HSIMT
      IF (Master) WRITE (stdout,20) 'TS_HSIMT',                        &
     &   'the 3rd HSIMT-TVD scheme for tracers (Wu and Zhu, 2010).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TS_HSIMT,'
      itrcHadv=itrcHadv+1
      itrcVadv=itrcVadv+1
# endif
# if defined TS_MIX_CLIMA && \
    (defined TS_DIF2      || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_MIX_CLIMA',                     &
     &   'Horizontal diffusion of tracer perturbation (T-Tclm).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TS_MIX_CLIMA,'
# endif
# if defined TS_MIX_MAX_SLOPE && defined MIX_ISO_TS && \
    (defined TS_DIF2          || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_MIX_MAX_SLOPE',                 &
     &   'Maximum density slope in isopycnal diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' TS_MIX_MAX_SLOPE,'
# endif
# if defined TS_MIX_MAX_SLOPE && defined MIX_ISO_TS && \
    (defined TS_DIF2          || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_MIX_MAX_STRAT',                 &
     &   'Minimum density stratification in isopycnal diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' TS_MIX_MAX_STRAT,'
# endif
# if defined TS_MIX_STABILITY && \
    (defined TS_DIF2          || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_MIX_STABILITY',                 &
     &   'Weighting diffusion between two time levels.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' TS_MIX_STABILITY,'
# endif
# ifdef TS_MPDATA
      IF (Master) WRITE (stdout,20) 'TS_MPDATA',                        &
     &   'Recursive flux corrected MPDATA 3D advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' TS_MPDATA,'
      itrcHadv=itrcHadv+1
      itrcVadv=itrcVadv+1
#  ifdef TS_MPDATA_LIMIT
      IF (Master) WRITE (stdout,20) 'TS_MPDATA_LIMIT',                  &
     &   'Further limiter in upwind corrector fluxes for stability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_MPDATA_LIMIT,'
#  endif
# endif
# ifdef TS_A4HADVECTION
      IF (Master) WRITE (stdout,20) 'TS_A4HADVECTION',                  &
     &   'Fourth-order Akima horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_A4HADVECTION,'
      itrcHadv=itrcHadv+1
# endif
# if defined TS_A4HADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_A4HADVECTION_TL',               &
     &   'TL/AD fourth-order Akima horizontal tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_A4HADVECTION_TL,'
      itrcHadvtl=itrcHadvtl+1
# endif
# ifdef TS_C2HADVECTION
      IF (Master) WRITE (stdout,20) 'TS_C2HADVECTION',                  &
     &   'Second-order centered horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C2HADVECTION,'
      itrcHadv=itrcHadv+1
# endif
# if defined TS_C2HADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_C2HADVECTION_TL',               &
     &   'TL/AD second-order centered horizontal tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_C2HADVECTION_TL,'
      itrcHadvtl=itrcHadvtl+1
# endif
# if defined TS_C4HADVECTION || defined TS_U3ADV_SPLIT
      IF (Master) WRITE (stdout,20) 'TS_C4HADVECTION',                  &
     &   'Fourth-order centered horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4HADVECTION,'
      itrcHadv=itrcHadv+1
# endif
# if defined TS_C4HADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_C4HADVECTION_TL',               &
     &   'TL/AD fourth-order centered horizontal tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_C4HADVECTION_TL,'
      itrcHadvtl=itrcHadvtl+1
# endif
# ifdef TS_U3HADVECTION
      IF (Master) WRITE (stdout,20) 'TS_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_U3HADVECTION,'
      itrcHadv=itrcHadv+1
# endif
# if defined TS_U3HADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_U3HADVECTION_TL',               &
     &   'TL/AD third-order upstream horizontal tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_U3HADVECTION_TL,'
      itrcHadvtl=itrcHadvtl+1
# endif
# if !(defined TS_MPDATA       || defined TS_A4HADVECTION || \
       defined TS_C2HADVECTION || defined TS_C4HADVECTION || \
       defined TS_U3HADVECTION || defined TS_HSIMT) 
      IF (Master) WRITE (stdout,20) 'TS_C4HADVECTION',                  &
     &   'Fourth-order centered horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4HADVECTION,'
      itrcHadv=itrcHadv+1
# endif
# ifdef TS_A4VADVECTION
      IF (Master) WRITE (stdout,20) 'TS_A4VADVECTION',                  &
     &   'Fourth-order Akima vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_A4VADVECTION,'
      itrcVadv=itrcVadv+1
# endif
# if defined TS_A4VADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_A4VADVECTION_TL',               &
     &   'TL/AD fourth-order Akima vertical tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_A4VADVECTION_TL,'
      itrcVadvtl=itrcVadvtl+1
# endif
# ifdef TS_C2VADVECTION
      IF (Master) WRITE (stdout,20) 'TS_C2VADVECTION',                  &
     &   'Second-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C2VADVECTION,'
      itrcVadv=itrcVadv+1
# endif
# if defined TS_C2VADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_C2VADVECTION_TL',               &
     &   'TL/AD second-order centered vertical tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_C2VADVECTION_TL,'
      itrcVadvtl=itrcVadvtl+1
# endif
# if defined TS_C4VADVECTION || defined TS_U3ADV_SPLIT
      IF (Master) WRITE (stdout,20) 'TS_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4VADVECTION,'
      itrcVadv=itrcVadv+1
# endif
# if defined TS_C4VADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_C4VADVECTION_TL',               &
     &   'TL/AD fourth-order centered vertical tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' TS_C4VADVECTION_TL,'
      itrcVadvtl=itrcVadvtl+1
# endif
# ifdef TS_SVADVECTION
      IF (Master) WRITE (stdout,20) 'TS_SVADVECTION',                   &
     &   'Parabolic splines vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' TS_SVADVECTION,'
      itrcVadv=itrcVadv+1
# endif
# if defined TS_SVADVECTION_TL && \
    (defined TANGENT || defined ADJOINT || defined TL_IOMS)
      IF (Master) WRITE (stdout,20) 'TS_SVADVECTION_TL',                &
     &   'TL/AD parabolic splines vertical tracer advection.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' TS_SVADVECTION_TL,'
      itrcVadvtl=itrcVadvtl+1
# endif
# if !(defined TS_MPDATA       || defined TS_A4VADVECTION || \
       defined TS_C2VADVECTION || defined TS_C4VADVECTION || \
       defined TS_SVADVECTION  || defined TS_HSIMT)
      IF (Master) WRITE (stdout,20) 'TS_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4VADVECTION,'
      itrcVadv=itrcVadv+1
# endif
# ifdef TS_U3ADV_SPLIT
      IF (Master) WRITE (stdout,20) 'TS_U3ADV_SPLIT',                   &
     &   'Split third-order upstream advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' TS_U3ADV_SPLIT,'
# endif
#endif
#if defined TS_DIF2 && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TS_DIF2',                          &
     &   'Harmonic mixing of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF2,'
#endif
#if defined TS_DIF4 && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TS_DIF4',                          &
     &   'Biharmonic mixing of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF4,'
#endif
#ifdef TS_FIXED
      IF (Master) WRITE (stdout,20) 'TS_FIXED',                         &
     &   'Diagnostic configuration, no evolution of tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TS_FIXED,'
#endif
#if defined TS_SMAGORINSKY && (defined TS_DIF2 || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_SMAGORINSKY',                   &
     &   'Smagorinksy-like time-dependent diffusion coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' TS_SMAGORINSKY,'
#endif
#if defined TWO_D_TRACER_SOURCE
      IF (Master) WRITE (stdout,20) 'TWO_D_TRACER_SOURCE',              &
     &   'One value per tracer per source.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' TWO_D_TRACER_SOURCE,'
#endif
#ifdef UV_ADV
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
#endif
#ifdef UV_COR
      IF (Master) WRITE (stdout,20) 'UV_COR',                           &
     &   'Coriolis term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
#endif
#ifdef UV_ADV
# ifdef UV_C2ADVECTION
      IF (Master) WRITE (stdout,20) 'UV_C2ADVECTION',                   &
     &   'Second-order centered differences advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C2ADVECTION,'
      ivelHadv=ivelHadv+1
      ivelVadv=ivelVadv+1
# endif
# if defined UV_C4ADVECTION || defined UV_U3ADV_SPLIT
      IF (Master) WRITE (stdout,20) 'UV_C4ADVECTION',                   &
     &   'Fourth-order centered differences advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4ADVECTION,'
      ivelHadv=ivelHadv+1
      ivelVadv=ivelVadv+1
# endif
# ifdef SOLVE3D
#  if !(defined UV_C2ADVECTION  || defined UV_C4ADVECTION)
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
#  endif
#  if !(defined UV_C2ADVECTION  || defined UV_C4ADVECTION || \
        defined UV_SADVECTION)
      IF (Master) WRITE (stdout,20) 'UV_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4VADVECTION,'
      ivelVadv=ivelVadv+1
#  endif
#  ifdef UV_KIRBY
      IF (Master) WRITE (stdout,20) 'UV_KIRBY',                         &
     &   'Compute uwave and vwave Kirby avg velocities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_KIRBY,'
#  endif
#  ifdef UV_SADVECTION
      IF (Master) WRITE (stdout,20) 'UV_SADVECTION',                    &
     &   'Parabolic splines vertical advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_SADVECTION,'
      ivelVadv=ivelVadv+1
#  endif
#  ifdef UV_U3ADV_SPLIT
      IF (Master) WRITE (stdout,20) 'UV_U3ADV_SPLIT',                   &
     &   'Split third-order upstream advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_U3ADV_SPLIT,'
#  endif
# endif
#endif
#ifdef UV_DRAG_GRID
      IF (Master) WRITE (stdout,20) 'UV_DRAG_GRID',                     &
# if defined UV_LOGDRAG
     &   'Spatially varying bottom roughness length.'
# elif defined UV_LDRAG
     &   'Spatially varying linear drag coefficient.'
# elif defined UV_QDRAG
     &   'Spatially varying quadratic drag coefficient.'
# endif
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' UV_DRAG_GRID,'
#endif
#ifdef UV_LDRAG
      IF (Master) WRITE (stdout,20) 'UV_LDRAG',                         &
     &   'Linear bottom stress.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_LDRAG,'
      ibbl=ibbl+1
#endif
#ifdef UV_LOGDRAG
      IF (Master) WRITE (stdout,20) 'UV_LOGDRAG',                       &
     &   'Logarithmic bottom stress.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' UV_LOGDRAG,'
      ibbl=ibbl+1
#endif
#ifdef UV_QDRAG
      IF (Master) WRITE (stdout,20) 'UV_QDRAG',                         &
     &   'Quadratic bottom stress.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_QDRAG,'
      ibbl=ibbl+1
#endif
#if defined UV_TIDES
      IF (Master) WRITE (stdout,20) 'UV_TIDES',                         &
     &   'Add tidal currents to 2D momentum climatologies.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_TIDES,'
#endif
#ifdef UV_VIS2
      IF (Master) WRITE (stdout,20) 'UV_VIS2',                          &
     &   'Harmonic mixing of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS2,'
#endif
#ifdef UV_VIS4
      IF (Master) WRITE (stdout,20) 'UV_VIS4',                          &
     &   'Biharmonic mixing of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS4,'
#endif
#if defined UV_SMAGORINSKY && (defined UV_VIS2 || defined UV_VIS4)
      IF (Master) WRITE (stdout,20) 'UV_SMAGORINSKY',                   &
     &   'Smagorinksy-like time-dependent viscosity coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_SMAGORINSKY,'
#endif
#ifdef UV_WAVEDRAG
      IF (Master) WRITE (stdout,20) 'UV_WAVEDRAG',                      &
     &   'Linear bottom stress from internal tidal wave drag.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' UV_WAVEDRAG,'
#endif
#if defined VAR_RHO_2D && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
#endif
#if defined UV_VIS2 || defined UV_VIS4
# ifdef VISC_GRID
      IF (Master) WRITE (stdout,20) 'VISC_GRID',                        &
     &   'Horizontal viscosity coefficient scaled by grid size.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' VISC_GRID,'
# endif
#endif
#if defined VCONVOLUTION && defined FOUR_DVAR && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'VCONVOLUTION',                     &
     &   'Include vertical correlations in convolutions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' VCONVOLUTION,'
#endif
#ifdef VERIFICATION
      IF (Master) WRITE (stdout,20) 'VERIFICATION',                     &
     &   'Proccess model solution at observation locations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' VERIFICATION,'
#endif
#if defined FLOATS && defined FLOAT_VWALK && defined SOLVE3D
# ifdef VWALK_FORWARD
      IF (Master) WRITE (stdout,20) 'VWALK_FORWARD',                    &
     &   'Lagrangian drifters, forward stepping random walk.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' VWALK_FORWARD,'
# else
      IF (Master) WRITE (stdout,20) '!VWALK_FORWARD',                   &
     &   'Lagrangian drifters, predictor/corrector random walk.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' !VWALK_FORWARD,'
# endif
#endif
#if defined VISC_3DCOEF
      IF (Master) WRITE (stdout,20) 'VISC_3DCOEF',                      &
     &   'Horizontal, time-dependent 3D viscosity coefficient.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' VISC_3DCOEF,'
#endif
#ifdef W4DPSAS
      IF (Master) WRITE (stdout,20) 'W4DPSAS',                          &
     &   'Weak constraint 4D-PSAS data assimilation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' W4DPSAS,'
      idriver=idriver+1
#endif
#ifdef W4DPSAS_SENSITIVITY
      IF (Master) WRITE (stdout,20) 'W4DPSAS_SENSITIVITY',              &
     &   'Weak constraint 4D-PSAS observation sensitivity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' W4DPSAS_SENSITIVITY,'
      idriver=idriver+1
#endif
#ifdef W4DVAR
      IF (Master) WRITE (stdout,20) 'W4DVAR',                           &
     &   'Weak constraint 4DVAR data assimilation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' W4DVAR,'
      idriver=idriver+1
#endif
#ifdef W4DVAR_SENSITIVITY
      IF (Master) WRITE (stdout,20) 'W4DVAR_SENSITIVITY',               &
     &   'Weak constraint 4DVAR observation sensitivity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' W4DVAR_SENSITIVITY,'
      idriver=idriver+1
#endif
#ifdef WAVES_OCEAN
      IF (Master) WRITE (stdout,20) 'WAVES_OCEAN',                      &
     &   'Two-way wave-ocean models coupling.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' WAVES_OCEAN,'
#endif
#ifdef WEAK_CONSTRAINT
      IF (Master) WRITE (stdout,20) 'WEAK_CONSTRAINT',                  &
     &   'Activated weak constraint assimilation set-up.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEAK_CONSTRAINT,'
#endif
#if defined WRITE_WATER && defined MASKING
      IF (Master) WRITE (stdout,20) 'WRITE_WATER',                      &
     &   'Writing data at water points only.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' WRITE_WATER,'
#endif
#ifdef WET_DRY
      IF (Master) WRITE (stdout,20) 'WET_DRY',                          &
     &   'Wetting and drying activated.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' WET_DRY,'
#endif
#ifdef WIND_WAVES
      IF (Master) WRITE (stdout,20) 'WIND_WAVES',                       &
     &   'Wind-induced wave data available.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' WIND_WAVES,'
#endif
#if !(defined PJ_GRADPQ4 || defined PJ_GRADPQ2 || defined PJ_GRADP || \
      defined DJ_GRADPS) && defined WJ_GRADP
      IF (Master) WRITE (stdout,20) 'WJ_GRADP',                         &
     &   'Weighted Jacobians pressure gradient adjustment.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' WJ_GRADP,'
#endif
#ifdef WRF_COUPLING
      IF (Master) WRITE (stdout,20) 'WRF_COUPLING',                     &
     &   'WRF model coupling.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' WRF_COUPLING,'
#endif
#if defined WTYPE_GRID && \
   (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
   !defined ANA_WTYPE
      IF (Master) WRITE (stdout,20) 'WTYPE_GRID',                       &
     &   'Spatially varying Jerlov water type index.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' WTYPE_GRID,'
#endif
#if defined ZETA_ELLIPTIC && defined BALANCE_OPERATOR && \
    defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'ZETA_ELLIPTIC',                    &
     &   'Solving SSH elliptic equation in balance operator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ZETA_ELLIPTIC,'
#endif
#ifdef ZOS_HSIG
      IF (Master) WRITE (stdout,20) 'ZOS_HSIG',                         &
     &   'Wave amplitude used for Zos calculation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ZOS_HSIG,'
#endif
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF

#ifndef SOLVE3D
!
!  Stop it explicit time-step splitting on shallow water set-up.
!
      DO ng=1,Ngrids
        IF (Master.and.(ndtfast(ng).gt.1)) THEN
          WRITE (stdout,80)
  80      FORMAT (/,' CHECKDEFS - explicit time-step splitting is ',    &
     &              ' inconsistent.',                                   &
     &            /,13x,'Change parameter NDTFAST to unity.')
          exit_flag=5
        END IF
      END DO
#endif
#if defined ADJOINT && defined _OPENMP
      IF (Master) THEN
        WRITE (stdout,90)
  90    FORMAT (/,' CHECKDEFS - cannot run the adjoint model in',       &
     &            ' shared-memory mode.')
        exit_flag=5
      END IF
#endif
#if defined NEMURO && defined HOLLING_GRAZING && defined IVLEV_EXPLICIT
!
!  Stop if using more than one grazing algorithm.
!
      IF (Master) THEN
        WRITE (stdout,100)
 100    FORMAT (/,' CHECKDEFS - More than one grazing algorithm',       &
     &            ' selected for Nemuro model.')
        exit_flag=5
      END IF
#endif
#if defined FOUR_DVAR && defined VCONVOLUTION
# if defined IMPLICIT_VCONV && defined SPLINES_VCONV
!
!  Stop if using more than vertical convolution algorithm.
!
      IF (Master) THEN
        WRITE (stdout,110)
 110    FORMAT (/,' CHECKDEFS - More than one vertical convolution',    &
     &            ' algorithm selected.')
        exit_flag=5
      END IF
# endif
#endif
#if defined COAWST_COUPLING && (defined ESMF_LIB && defined MCT_LIB)
      IF (Master) THEN
        WRITE (stdout,130)
 130    FORMAT (/,' CHECKDEFS - More than one coupling library',        &
     &            ' selected.')
        exit_flag=5
      END IF
#endif
#if defined BULK_FLUXES && defined ATM2OCN_FLUXES
      IF (Master) THEN
        WRITE (stdout,135)
 135    FORMAT (/,' CHECKDEFS - More than one surface flux method ',    &
     &            ' selected.')
        exit_flag=5
      END IF
#endif
#if defined POSTERIOR_ERROR_F && defined POSTERIOR_ERROR_I && \
    defined WEAK_CONSTRAINT
      IF (Master) THEN
        WRITE (stdout,130)
 130    FORMAT (/,' CHECKDEFS - Both initial and final posterior',      &
     &            ' error covariance analysis selected, choose one.')
        exit_flag=5
      END IF
#endif
#ifdef UV_ADV
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
# ifdef SOLVE3D
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
# endif
#endif
#ifdef SOLVE3D
      IF (Master.and.(itrcHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','tracers','itrcHadv =',itrcHadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','tracers','itrcVadv =',itrcVadv
        exit_flag=5
      END IF
# if defined TANGENT || defined TL_IOMS || defined ADJOINT
      IF (Master.and.(itrcHadvtl.gt.1)) THEN
        WRITE (stdout,140) 'TL/AD horizontal','tracers','itrcHadvtl =', &
     &                     itrcHadvtl
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadvtl.gt.1)) THEN
        WRITE (stdout,140) 'TL/AD vertical','tracers','itrcVadvtl =',   &
     &                     itrcVadvtl
        exit_flag=5
      END IF
# endif
#endif
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
#if defined SOLVE3D && defined NEMURO
# if !(defined IVLEV_EXPLICIT || defined HOLLING_GRAZING)
      IF (Master) THEN
        WRITE (stdout,150) uppercase('ivlev_explicit'),                 &
     &                     uppercase('holling_grazing')
 150    FORMAT (/,' CHECKDEFS - Need to choose a zooplankton grazing',  &
     &            ' option;',/,13x,'use either ',a,' or implicit ',a,   &
     &            '.',/,13x,'The default implicit IVLEV algorithm',     &
     &            ' does not work well yet.')
        exit_flag=5
      END IF
# endif
#endif
#ifdef SEDIMENT
      IF (Master.and.(sbed_type.gt.1)) THEN
        WRITE (stdout,170)
 170    FORMAT (/,' CHECKDEFS - More than one sediment bed type   ',    &
     &            ' selected.')
        exit_flag=5
      END IF
      IF (Master.and.(sbed_load.gt.1)) THEN
        WRITE (stdout,180)
 180    FORMAT (/,' CHECKDEFS - More than one sediment bed load   ',    &
     &            ' selected.')
        exit_flag=5
      END IF
#endif
#if defined ICESHELF && !defined SOLVE3D
!
!  Stop if ICESHELF defined in shallow water setup
!
      IF (Master) THEN
        WRITE (stdout,170)
 170    FORMAT (/,' CHECKDEFS - ICESHELF does not currently work',      &
     &            ' in shallow water mode.')
        exit_flag=5
      END IF
#endif
#if defined ICESHELF && !defined DJ_GRADPS
!
!  PGF for ICESHELF only works for DJ_GRADPS
!
      IF (Master) THEN
        WRITE (stdout,180)
 180    FORMAT (/,' CHECKDEFS - ICESHELF does not currently work',      &
     &            ' for any PGF scheme other than DJ_GRADPS')
#endif
#ifdef WEC
!
!  Stop if more that one radiation stress formulation is activated.
!
      IF (Master.and.(nearshore.gt.1)) THEN
        WRITE (stdout,160)
 160    FORMAT (/,' CHECKDEFS - only one wave current formulation'  &
     &            ' is allowed.')
        exit_flag=5
      END IF
#endif
#if defined PARALLEL_IO || defined PARALLEL_OUT
# if !(defined HDF5 || defined PNETCDF)
      IF (Master) THEN
        WRITE (stdout,170) uppercase('hdf5'),                           &
     &                     uppercase('pnetcdf')
 170    FORMAT (/,' CHECKDEFS - Need to activate either ',a,' or ',a,   &
     &            ' options',/,13x,'for parallel I/O processing.')
        exit_flag=5
      END IF
# endif
# if defined HDF5 && defined PNETCDF
      IF (Master) THEN
        WRITE (stdout,180) uppercase('hdf5'),                           &
     &                     uppercase('pnetcdf')
 180    FORMAT (/,' CHECKDEFS - Both ',a,' and ',a,' options',          &
     &            ' selected',/,13x,'for parallel I/O processing.',     &
     &          /,13x,'Choose only one.')
        exit_flag=5
      END IF
# endif
#endif
#if defined WRITE_WATER && defined WET_DRY
!
!  Stop if writting only water points when wetting and drying is
!  activated.  This is not possible since the Land/Sea masking is
!  time dependent.
!
      IF (Master) THEN
        WRITE (stdout,190) uppercase('write_water'),                    &
     &                     uppercase('wet_dry')
 190    FORMAT (/,' CHECKDEFS - cannot activate ',a,' and ',a,          &
     &            ' together',/,13x,'because of time dependent ',       &
     &            ' Land/Sea masking.')
        exit_flag=5
      END IF
#endif
#if ((defined AD_AVERAGES && defined ADJOINT) && defined AVERAGES) || \
    ((defined RP_AVERAGES && defined TL_IOMS) && defined AVERAGES) || \
    ((defined TL_AVERAGES && defined TANGENT) && defined AVERAGES) || \
    (defined AD_AVERAGES && defined TL_AVERAGES) || \
    (defined AD_AVERAGES && defined RP_AVERAGES) || \
    (defined RP_AVERAGES && defined TL_AVERAGES)
!
!  Stop if activating more that one time-averaged option. This is not
!  possible in the same executable because the same internal arrays
!  are used to store the cummulative sums.
!
      IF (Master) THEN
        WRITE (stdout,200) uppercase('averages'),                       &
     &                     uppercase('ad_averages'),                    &
     &                     uppercase('rp_averages'),                    &
     &                     uppercase('tl_averages')
 200    FORMAT (/,' CHECKDEFS - cannot activate ',a,' ',a,' ',a,        &
     &            ' or',a,/,13x,'at the same time. Choose one!')
        exit_flag=5
      END IF
#endif
#if !defined DISTRIBUTE && defined TS_MPDATA
!
!  Stop if activating MPDATA in serial with partitions or shared-
!  memory.
!
      IF (Master) THEN
        DO ng=1,Ngrids
          IF (NtileX(ng)*NtileE(ng).gt.1) THEN
            WRITE (stdout,210) uppercase('ts_mpdata')
            exit_flag=5
            EXIT
          END IF
        END DO
 210    FORMAT (/,' CHECKDEFS - cannot activate option: ',a,           &
     &          /,13x,'in serial with partitions or shared-memory...', &
     &          /,13x,'Use distributed-memory (MPI) in parallel runs.')
      END IF
#endif
#ifdef Q_PSOURCE
!
!  Stop if activating obsolete mass point Sources/Sinks option. This
!  capability is now activated with standard input switch "LwSrc".
!
      IF (Master) THEN
        WRITE (stdout,220) uppercase('q_psource'), 'LwSrc'
 220    FORMAT (/,' CHECKDEFS - cannot use obsolete option: ',a,        &
     &          /,13x,'Use instead standard input switch: ',a,          &
     &          /,13x,'Edit header file or build script!')
        exit_flag=5
      END IF
#endif
#if defined TS_PSOURCE && defined SOLVE3D
!
!  Stop if activating obsolete tracers point Sources/Sinks option. This
!  capability is now activated with standard input switch "LtracerSrc".
!
      IF (Master) THEN
        WRITE (stdout,230) uppercase('ts_psource'), 'LtracerSrc'
 230    FORMAT (/,' CHECKDEFS - cannot use obsolete option: ',a,        &
     &          /,13x,'Use instead standard input switch: ',a,          &
     &          /,13x,'Edit header file or build script!')
        exit_flag=5
      END IF
#endif
#ifdef UV_PSOURCE
!
!  Stop if activating obsolete mass point Sources/Sinks option. This
!  capability is now activated with standard input switch "LwSrc".
!
      IF (Master) THEN
        WRITE (stdout,240) uppercase('uv_psource'), 'LuvSrc'
 240    FORMAT (/,' CHECKDEFS - cannot use obsolete option: ',a,        &
     &          /,13x,'Use instead standard input switch: ',a,          &
     &          /,13x,'Edit header file or build script!')
        exit_flag=5
      END IF
#endif
#if defined SPLINES && defined SOLVE3D
!
!  Stop if activating parabolic splines reconstruction with deprecated
!  option. See https://www.myroms.org/projects/src/ticket/681 for more
!  details.
!
      IF (Master) THEN
        WRITE (stdout,250) uppercase('splines'),                        &
     &        'https://www.myroms.org/projects/src/ticket/681'
 250    FORMAT (/,' CHECKDEFS - cannot use obsolete option: ',a,        &
     &          /,13x,'See following trac ticket for details: ',        &
     &          /,13x,a,/,13x,'please change header file accordingly.')
        exit_flag=5
      END IF
#endif
#if !defined DISTRIBUTE && defined NESTING
!
!  Stop if activating nesting in serial with partitions or shared-
!  memory.
!
      IF (Master) THEN
        DO ng=1,Ngrids
          IF (NtileX(ng)*NtileE(ng).gt.1) THEN
            WRITE (stdout,260) uppercase('nesting')
            exit_flag=5
            EXIT
          END IF
        END DO
 260    FORMAT (/,' CHECKDEFS - cannot activate option: ',a,           &
     &          /,13x,'in serial with partitions or shared-memory...', &
     &          /,13x,'There bugs in shared-memory that need fixing.', &
     &          /,13x,'Use distributed-memory (MPI) in parallel runs.')
      END IF
#endif
#if defined CICE_MODEL && defined ICE_MODEL
!
!  Stop if activating two sea ice models.
!
      IF (Master) THEN
        WRITE (stdout,260) uppercase('cice_model'),                     &
     &                     uppercase('ice_model')
 260    FORMAT (/,' CHECKDEFS - cannot activate ',a,'                   &
     &            ' and',a,/,13x,'at the same time. Choose one!')
        exit_flag=5
      END IF
#endif


      RETURN
      END SUBROUTINE checkdefs
