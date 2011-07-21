#include "cppdefs.h"
      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 365 2007-07-08 16:15:48Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      implicit none
!
!  Local variable declarations.
!
      integer :: FS_east = 0
      integer :: FS_west = 0
      integer :: FS_north = 0
      integer :: FS_south = 0
      integer :: M2_east = 0
      integer :: M2_west = 0
      integer :: M2_north = 0
      integer :: M2_south = 0
      integer :: M3_east = 0
      integer :: M3_west = 0
      integer :: M3_north = 0
      integer :: M3_south = 0
      integer :: M2s_east = 0
      integer :: M2s_west = 0
      integer :: M2s_north = 0
      integer :: M2s_south = 0
      integer :: M3s_east = 0
      integer :: M3s_west = 0
      integer :: M3s_north = 0
      integer :: M3s_south = 0
      integer :: K_east = 0
      integer :: K_west = 0
      integer :: K_north = 0
      integer :: K_south = 0
      integer :: T_east = 0
      integer :: T_west = 0
      integer :: T_north = 0
      integer :: T_south = 0
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
#if defined ADD_FSOBC && defined SSH_TIDES
      IF (Master) WRITE (stdout,20) 'ADD_FSOBC',                        &
     &   'Adding tidal elevation to proccesed OBC data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_FSOBC,'
#endif
#if defined ADD_M2OBC && defined UV_TIDES
      IF (Master) WRITE (stdout,20) 'ADD_M2OBC',                        &
     &   'Adding tidal currents to proccesed OBC data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_M2OBC,'
#endif
#ifdef ADJOINT
      IF (Master) WRITE (stdout,20) 'ADJOINT',                          &
     &   'Adjoint Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ADJOINT,'
#endif
#if defined ADJUST_BOUNDARY   && \
   (defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
    defined SENSITIVITY_4DVAR || defined W4DPSAS             || \
    defined W4DVAR)
      IF (Master) WRITE (stdout,20) 'ADJUST_BOUNDARY',                  &
     &   'Including boundary conditions in 4DVar state estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ADJUST_BOUNDARY,'
#endif
#if defined ADJUST_STFLUX     && \
   (defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
    defined SENSITIVITY_4DVAR || defined W4DPSAS             || \
    defined W4DVAR)
      IF (Master) WRITE (stdout,20) 'ADJUST_STFLUX',                    &
     &   'Including surface tracer flux in 4DVar state estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ADJUST_STFLUX,'
#endif
#if defined ADJUST_WSTRESS    && \
   (defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
    defined SENSITIVITY_4DVAR || defined W4DPSAS             || \
    defined W4DVAR)
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
#if defined ALBEDO && defined ANA_SRFLUX
      IF (Master) WRITE (stdout,20) 'ALBEDO',                           &
     &   'Shortwave radiation from albedo equation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ALBEDO,'
#endif
#ifdef ATM_PRESS
      IF (Master) WRITE (stdout,20) 'ATM_PRESS',                        &
     &   'Impose atmospheric pressure onto sea surface.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ATM_PRESS,'
#endif
#ifdef BALANCE_OPERATOR
      IF (Master) WRITE (stdout,20) 'BALANCE_OPERATOR',                 &
     &   'Error Covariance Multivariate Balance Operator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' BALANCE_OPERATOR,'
#endif
#ifdef CONVOLUTION
      IF (Master) WRITE (stdout,20) 'CONVOLUTION',                      &
     &   'Generic Convolution of State Solution Driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' CONVOLUTION,'
      idriver=idriver+1
#endif
#if defined CONVOLVE && defined W4DVAR && defined W4DVAR_SENSITIVITY
      IF (Master) WRITE (stdout,20) 'CONVOLVE',                         &
     &   'Convolve Adjoint State Solution.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CONVOLVE,'
#endif
#ifdef CORRELATION
      IF (Master) WRITE (stdout,20) 'CORRELATION',                      &
     &   'Background-error Correlation Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' CORRELATION,'
      idriver=idriver+1
#endif
#if defined BIOLOGY && defined ANA_BIOLOGY
      IF (Master) WRITE (stdout,20) 'ANA_BIOLOGY',                      &
     &   'Analytical biology initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_BIOLOGY,'
#endif
#if defined SEDIMENT && !defined ANA_SEDIMENT
# ifdef ANA_BMFLUX
      IF (Master) WRITE (stdout,20) 'ANA_BMFLUX',                       &
     &   'Analytical spatially varying bottom roughness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BMFLUX,'
# endif
#endif
#if defined SEDIMENT || defined BIOLOGY
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
#ifdef ANA_DIAG
      IF (Master) WRITE (stdout,20) 'ANA_DIAG',                         &
     &   'Customized diagnostics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_DIAG,'
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
#ifdef ANA_INITIAL
      IF (Master) WRITE (stdout,20) 'ANA_INITIAL',                      &
     &   'Analytical initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_INITIAL,'
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
#ifdef ANA_SMFLUX
      IF (Master) WRITE (stdout,20) 'ANA_SMFLUX',                       &
     &   'Analytical kinematic surface momentum flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SMFLUX,'
#endif
#if defined SEDIMENT || defined BIOLOGY
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
#ifdef ANA_SST
      IF (Master) WRITE (stdout,20) 'ANA_SST',                          &
     &   'Analytical SST and dQdSST.'
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
#ifdef ANA_WWAVE
      IF (Master) WRITE (stdout,20) 'ANA_WWAVE',                        &
     &   'Analytical wind induced waves.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WWAVE,'
#endif
#ifdef ASSIMILATION_SSH
      IF (Master) WRITE (stdout,20) 'ASSIMILATION_SSH',                 &
     &   'Assimilating SSH observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' ASSIMILATION_SSH,'
#endif
#ifdef ASSIMILATION_SST
      IF (Master) WRITE (stdout,20) 'ASSIMILATION_SST',                 &
     &   'Assimilating SST observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' ASSIMILATION_SST,'
#endif
#ifdef ASSIMILATION_T
      IF (Master) WRITE (stdout,20) 'ASSIMILATION_T',                   &
     &   'Assimilating tracers observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ASSIMILATION_T,'
#endif
#ifdef ASSIMILATION_UV
      IF (Master) WRITE (stdout,20) 'ASSIMILATION_UV',                  &
     &   'Assimilating horizontal currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ASSIMILATION_UV,'
#endif
#ifdef ASSIMILATION_UVsur
      IF (Master) WRITE (stdout,20) 'ASSIMILATION_UVsur',               &
     &   'Assimilating surface currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' ASSIMILATION_UVsur,'
#endif
#ifdef ASSUMED_SHAPE
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
#endif
#ifdef ATM2OCN_FLUXES
      IF (Master) WRITE (stdout,20) 'ATM2OCN_FLUXES',                   &
     &   'Surface fluxes from ATM Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ATM2OCN_FLUXES,'
#endif
#ifdef AVERAGES
      IF (Master) WRITE (stdout,20) 'AVERAGES',                         &
     &   'Writing out time-averaged fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AVERAGES,'
# ifdef AVERAGES_AKS
      IF (Master) WRITE (stdout,20) 'AVERAGES_AKS',                     &
     &   'Writing out time-averaged vertical S-diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' AVERAGES_AKS,'
# endif
# ifdef AVERAGES_AKT
      IF (Master) WRITE (stdout,20) 'AVERAGES_AKT',                     &
     &   'Writing out time-averaged vertical T-diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' AVERAGES_AKT,'
# endif
# if defined AVERAGES_BEDLOAD && defined SEDIMENT && defined BEDLOAD
      IF (Master) WRITE (stdout,20) 'AVERAGES_BELOAD',                  &
     &   'Writing out time-averaged bed load.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' AVERAGES_BEDLOAD,'
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      IF (Master) WRITE (stdout,20) 'AVERAGES_DETIDE',                  &
     &   'Writing out time-averaged detided fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' AVERAGES_DETIDE,'
# endif
# ifdef AVERAGES_FLUXES
      IF (Master) WRITE (stdout,20) 'AVERAGES_FLUXES',                  &
     &   'Writing out time-averaged surface fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' AVERAGES_FLUXES,'
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
# ifdef AVERAGES_QUADRATIC
      IF (Master) WRITE (stdout,20) 'AVERAGES_QUADRATIC',               &
     &   'Writing out time-averaged quadratic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' AVERAGES_QUADRATIC,'
# endif
#endif
#if defined BACKGROUND && defined FOUR_DVAR
      IF (Master) WRITE (stdout,20) 'BACKGROUND',                       &
     &   'Include bacground cost function.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BACKGROUND,'
#endif
#if defined SEDIMENT && defined BEDLOAD_MPM
      IF (Master) WRITE (stdout,20) 'BEDLOAD_MPM',                      &
     &   'Activate bed load sediment transport Meyer-Peter-Mueller.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BEDLOAD_MPM,'
#endif
#if defined SEDIMENT && defined BEDLOAD_SOULSBY
      IF (Master) WRITE (stdout,20) 'BEDLOAD_SOULSBY',                  &
     &   'Activate bed load sediment transport Soulsby formula.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' BEDLOAD_SOULSBY,'
#endif
#ifdef BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'BIO_FENNEL',                       &
     &   'Fennel et al. (2006) nitrogen-based model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BIO_FENNEL,'
      ibiology=ibiology+1
#endif
#if defined BIO_SEDIMENT && \
   (defined BIO_FENNEL   || defined NEMURO || defined ECOSIM)
      IF (Master) WRITE (stdout,20) 'BIO_SEDIMENT',                     &
     &   'Restore fallen particulate material to the nutrient pool.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' BIO_SEDIMENT,'
#endif
#ifdef BODYFORCE
      IF (Master) WRITE (stdout,20) 'BODYFORCE',                        &
     &   'Momentum stresses as body-forces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' BODYFORCE,'
#endif
#ifdef BULK_FLUXES
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES',                      &
     &   'Surface bulk fluxes parametererization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
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
      Coptions(is:is+16)=' CELERETY_WRITE,'
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
#if defined DEFLATE && defined NETCDF4
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
#if defined DIAGNOSTICS_BIO && defined BIO_FENNEL
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_BIO',                  &
     &   'Computing and writing biological diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DIAGNOSTICS_BIO,'
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
#ifdef EAST_FSCLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_FSCLAMPED',                   &
     &   'Eastern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSCLAMPED,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSCHAPMAN
      IF (Master) WRITE (stdout,20) 'EAST_FSCHAPMAN',                   &
     &   'Eastern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSCHAPMAN,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSGRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_FSGRADIENT',                  &
     &   'Eastern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_FSGRADIENT,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSRADIATION
# ifdef EAST_FSNUDGING
      IF (Master) WRITE (stdout,20) 'EAST_FSNUDGING',                   &
     &   'Eastern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'EAST_FSRADIATION',                 &
     &   'Eastern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_FSRADIATION,'
      FS_east=FS_east+1
#endif
#ifdef EAST_M2CLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_M2CLAMPED',                   &
     &   'Eastern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2CLAMPED,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2FLATHER
      IF (Master) WRITE (stdout,20) 'EAST_M2FLATHER',                   &
     &   'Eastern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2FLATHER,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2GRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_M2GRADIENT',                  &
     &   'Eastern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M2GRADIENT,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2RADIATION
# ifdef EAST_M2NUDGING
      IF (Master) WRITE (stdout,20) 'EAST_M2NUDGING',                   &
     &   'Eastern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'EAST_M2RADIATION',                 &
     &   'Eastern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M2RADIATION,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2REDUCED
      IF (Master) WRITE (stdout,20) 'EAST_M2REDUCED',                   &
     &   'Eastern edge, 2D momentum, reduced-physics condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2REDUCED,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M3CLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_M3CLAMPED',                   &
     &   'Eastern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M3CLAMPED,'
      M3_east=M3_east+1
#endif
#ifdef EAST_M3GRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_M3GRADIENT',                  &
     &   'Eastern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M3GRADIENT,'
      M3_east=M3_east+1
#endif
#ifdef EAST_M3RADIATION
# ifdef EAST_M3NUDGING
      IF (Master) WRITE (stdout,20) 'EAST_M3NUDGING',                   &
     &   'Eastern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M3NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'EAST_M3RADIATION',                 &
     &   'Eastern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M3RADIATION,'
      M3_east=M3_east+1
#endif
#ifdef EAST_M2SCLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_M2SCLAMPED',                  &
     &   'Eastern edge, 2D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M2SCLAMPED,'
      M2s_east=M2s_east+1
#endif
#ifdef EAST_M2SGRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_M2SGRADIENT',                 &
     &   'Eastern edge, 2D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M2SGRADIENT,'
      M2s_east=M2s_east+1
#endif
#ifdef EAST_M2SRADIATION
      IF (Master) WRITE (stdout,20) 'EAST_M2SRADIATION',                &
     &   'Eastern edge, 2D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' EAST_M2SRADIATION,'
      M2s_east=M2s_east+1
#endif
#ifdef EAST_M3SCLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_M3SCLAMPED',                  &
     &   'Eastern edge, 3D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M3SCLAMPED,'
      M3s_east=M3s_east+1
#endif
#ifdef EAST_M3SGRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_M3SGRADIENT',                 &
     &   'Eastern edge, 3D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M3SGRADIENT,'
      M3s_east=M3s_east+1
#endif
#ifdef EAST_M3SRADIATION
      IF (Master) WRITE (stdout,20) 'EAST_M3SRADIATION',                 &
     &   'Eastern edge, 3D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' EAST_M3SRADIATION,'
      M3s_east=M3s_east+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef EAST_KGRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_KGRADIENT',                   &
     &   'Eastern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_KGRADIENT,'
      K_east=K_east+1
# endif
# ifdef EAST_KRADIATION
      IF (Master) WRITE (stdout,20) 'EAST_KRADIATION',                  &
     &   'Eastern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_KRADIATION,'
      K_east=K_east+1
# endif
#endif
#ifdef EAST_TCLAMPED
      IF (Master) WRITE (stdout,20) 'EAST_TCLAMPED',                    &
     &   'Eastern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' EAST_TCLAMPED,'
      T_east=T_east+1
#endif
#ifdef EAST_TGRADIENT
      IF (Master) WRITE (stdout,20) 'EAST_TGRADIENT',                   &
     &   'Eastern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_TGRADIENT,'
      T_east=T_east+1
#endif
#ifdef EAST_TRADIATION
# ifdef EAST_TNUDGING
      IF (Master) WRITE (stdout,20) 'EAST_TNUDGING',                    &
     &   'Eastern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' EAST_TNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'EAST_TRADIATION',                  &
     &   'Eastern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_TRADIATION,'
      T_east=T_east+1
#endif
#ifdef EAST_VOLCONS
      IF (Master) WRITE (stdout,20) 'EAST_VOLCONS',                     &
     &   'Eastern edge, enforce mass conservation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' EAST_VOLCONS,'
#endif
#ifdef EASTERN_WALL
      IF (Master) WRITE (stdout,20) 'EASTERN_WALL',                     &
     &   'Wall boundary at Eastern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' EASTERN_WALL,'
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
#if defined ESMF_LIB && defined MODEL_COUPLING
      IF (Master) WRITE (stdout,20) 'ESMF_LIB',                         &
     &   'Using Earth System Modeling Framework library.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ESMF_LIB,'
#endif
#ifdef ENSEMBLE
      IF (Master) WRITE (stdout,20) 'ENSEMBLE',                         &
     &   'Ensemble Forecasting Propagator.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ENSEMBLE,'
#endif
#ifdef EW_PERIODIC
      IF (Master) WRITE (stdout,20) 'EW_PERIODIC',                      &
     &   'East-West periodic boundaries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' EW_PERIODIC,'
      FS_east=FS_east+1
      FS_west=FS_west+1
      M2_east=M2_east+1
      M2_west=M2_west+1
      M3_east=M3_east+1
      M3_west=M3_west+1
# ifdef WEC_VF
      M2s_east=M2s_east+1
      M2s_west=M2s_west+1
      M3s_east=M3s_east+1
      M3s_west=M3s_west+1
# endif
      K_east=K_east+1
      K_west=K_west+1
      T_east=T_east+1
      T_west=T_west+1
#endif
#ifdef FLOATS
      IF (Master) WRITE (stdout,20) 'FLOATS',                           &
     &   'Simulated Lagrangian drifters.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' FLOATS,'
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
     &   'Using free-surface data in reduced physics condtions'
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
#ifdef GLS_MIXING
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
#endif
#if defined HOLLING_GRAZING && defined NEMURO
      IF (Master) WRITE (stdout,20) 'HOLLING_GRAZING',                  &
     &   'Holling-type s-shaped grazing algorithm, Implicit.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' HOLLING_GRAZING,'
#endif
#ifdef ICESHELF
      IF (Master) WRITE (stdout,20) 'ICESHELF',                         &
     &   'Include Ice Shelf Cavities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ICESHELF,'
#endif
#if defined IMPLICIT_VCONV && defined FOUR_DVAR && defined VCONVOLUTION
      IF (Master) WRITE (stdout,20) 'IMPLICIT_VCONV',                   &
     &   'Implicit Vertical Convolution Algorithm.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' IMPLICIT_VCONV,'
#endif
#if defined NPZD_IRON && defined IRON_LIMIT
      IF (Master) WRITE (stdout,20) 'IRON_LIMIT',                       &
     &   'Include dissolved iron limitation on phytoplankton growth.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' IRON_LIMIT,'
#endif
#if defined NPZD_IRON && defined IRON_LIMIT && defined IRON_RELAX
      IF (Master) WRITE (stdout,20) 'IRON_RELAX',                       &
     &   'Nudge dissolved iron over the shelf, h <= FeHmin.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' IRON_RELAX,'
#endif
#ifdef IMPULSE
      IF (Master) WRITE (stdout,20) 'IMPULSE',                          &
     &   'Processing Adjoint Impulse forcing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' IMPULSE,'
#endif
#ifdef INNER_PRODUCT
      IF (Master) WRITE (stdout,20) 'INNER_PRODUCT',                    &
     &   'Tangent/Adjoint State Matrices Inner Product Test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' INNER_PRODUCT,'
      idriver=idriver+1
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
#ifdef M2CLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'M2CLIMATOLOGY',                    &
     &   'Processing 2D momentum climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M2CLIMATOLOGY,'
#endif
#ifdef M2CLM_NUDGING
      IF (Master) WRITE (stdout,20) 'M2CLM_NUDGING',                    &
     &   'Nudging toward 2D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M2CLM_NUDGING,'
#endif
#if defined M3CLIMATOLOGY && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'M3CLIMATOLOGY',                    &
     &   'Processing 3D momentum climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M3CLIMATOLOGY,'
#endif
#if defined M3CLM_NUDGING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'M3CLM_NUDGING',                    &
     &   'Nudging toward 3D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M3CLM_NUDGING,'
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
#if defined MCT_LIB && defined MODEL_COUPLING
      IF (Master) WRITE (stdout,20) 'MCT_LIB',                          &
     &   'Using Model Coupling Toolkit library.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MCT_LIB,'
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
#endif
#ifdef WEC_VF
      IF (Master) WRITE (stdout,20) 'WEC_VF',                           &
     &   'Vortex Force wave current interaction- UNDER DEVELOPMENT!!!!.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' WEC_VF,'
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
#ifdef NETCDF4
      IF (Master) WRITE (stdout,20) 'NETCDF4',                          &
     &   'Creating NetCDF-4/HDF5 format files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' NETCDF4,'
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
#ifdef NORTH_VOLCONS
      IF (Master) WRITE (stdout,20) 'NORTH_VOLCONS',                    &
     &   'Northern edge, enforce mass conservation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NORTH_VOLCONS,'
#endif
#ifdef NORTHERN_WALL
      IF (Master) WRITE (stdout,20) 'NORTHERN_WALL',                    &
     &   'Wall boundary at Northern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NORTHERN_WALL,'
#endif
#ifdef NORTH_FSCLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_FSCLAMPED',                  &
     &   'Northern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSCLAMPED,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSCHAPMAN
      IF (Master) WRITE (stdout,20) 'NORTH_FSCHAPMAN',                  &
     &   'Northern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSCHAPMAN,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSGRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_FSGRADIENT',                 &
     &   'Northern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_FSGRADIENT,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSRADIATION
# ifdef NORTH_FSNUDGING
      IF (Master) WRITE (stdout,20) 'NORTH_FSNUDGING',                  &
     &   'Northern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'NORTH_FSRADIATION',                &
     &   'Northern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_FSRADIATION,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_M2CLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_M2CLAMPED',                  &
     &   'Northern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2CLAMPED,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2FLATHER
      IF (Master) WRITE (stdout,20) 'NORTH_M2FLATHER',                  &
     &   'Northern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2FLATHER,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2GRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_M2GRADIENT',                 &
     &   'Northern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M2GRADIENT,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2RADIATION
# ifdef NORTH_M2NUDGING
      IF (Master) WRITE (stdout,20) 'NORTH_M2NUDGING',                  &
     &   'Northern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'NORTH_M2RADIATION',                &
     &   'Northern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M2RADIATION,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2REDUCED
      IF (Master) WRITE (stdout,20) 'NORTH_M2REDUCED',                  &
     &   'Northern edge, 2D momentum, reduced-physics condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2REDUCED,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M3CLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_M3CLAMPED',                  &
     &   'Northern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M3CLAMPED,'
      M3_north=M3_north+1
#endif
#ifdef NORTH_M3GRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_M3GRADIENT',                 &
     &   'Northern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M3GRADIENT,'
      M3_north=M3_north+1
#endif
#ifdef NORTH_M3RADIATION
# ifdef NORTH_M3NUDGING
      IF (Master) WRITE (stdout,20) 'NORTH_M3NUDGING',                  &
     &   'Northern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M3NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'NORTH_M3RADIATION',                &
     &   'Northern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M3RADIATION,'
      M3_north=M3_north+1
#endif
#ifdef NORTH_M2SCLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_M2SCLAMPED',                 &
     &   'Northern edge, 2D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M2SCLAMPED,'
      M2s_north=M2s_north+1
#endif
#ifdef NORTH_M2SGRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_M2SGRADIENT',                &
     &   'Northern edge, 2D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M2SGRADIENT,'
      M2s_north=M2s_north+1
#endif
#ifdef NORTH_M2SRADIATION
      IF (Master) WRITE (stdout,20) 'NORTH_M2SRADIATION',               &
     &   'Northern edge, 2D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' NORTH_M2SRADIATION,'
      M2s_north=M2s_north+1
#endif
#ifdef NORTH_M3SCLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_M3SCLAMPED',                 &
     &   'Northern edge, 3D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M3SCLAMPED,'
      M3s_north=M3s_north+1
#endif
#ifdef NORTH_M3SGRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_M3SGRADIENT',                &
     &   'Northern edge, 3D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M3SGRADIENT,'
      M3s_north=M3s_north+1
#endif
#ifdef NORTH_M3SRADIATION
      IF (Master) WRITE (stdout,20) 'NORTH_M3SRADIATION',                 &
     &   'Northern edge, 3D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' NORTH_M3SRADIATION,'
      M3s_north=M3s_north+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef NORTH_KGRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_KGRADIENT',                  &
     &   'Northern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_KGRADIENT,'
      K_north=K_north+1
# endif
# ifdef NORTH_KRADIATION
      IF (Master) WRITE (stdout,20) 'NORTH_KRADIATION',                 &
     &   'Northern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_KRADIATION,'
      K_north=K_north+1
# endif
#endif
#ifdef NORTH_TCLAMPED
      IF (Master) WRITE (stdout,20) 'NORTH_TCLAMPED',                   &
     &   'Northern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' NORTH_TCLAMPED,'
      T_north=T_north+1
#endif
#ifdef NORTH_TGRADIENT
      IF (Master) WRITE (stdout,20) 'NORTH_TGRADIENT',                  &
     &   'Northern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_TGRADIENT,'
      T_north=T_north+1
#endif
#ifdef NORTH_TRADIATION
# ifdef NORTH_TNUDGING
      IF (Master) WRITE (stdout,20) 'NORTH_TNUDGING',                   &
     &   'Northern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' NORTH_TNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'NORTH_TRADIATION',                 &
     &   'Northern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_TRADIATION,'
      T_north=T_north+1
#endif
#if defined NO_READ_GHOST && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'NO_READ_GRID',                     &
     &   'Not including ghost points during scattering read data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' NO_READ_GRID,'
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
     &   'NPZD Biological Model (Powell et al.) with iron limitation'
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
#ifdef NUDGING_SST
      IF (Master) WRITE (stdout,20) 'NUDGING_SST',                      &
     &   'Nudging SST observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NUDGING_SST,'
#endif
#ifdef NUDGING_T
      IF (Master) WRITE (stdout,20) 'NUDGING_T',                        &
     &   'Nudging tracers observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' NUDGING_T,'
#endif
#ifdef NUDGING_UV
      IF (Master) WRITE (stdout,20) 'NUDGING_UV',                       &
     &   'Nudging horizontal currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NUDGING_UV,'
#endif
#ifdef NUDGING_UVsur
      IF (Master) WRITE (stdout,20) 'NUDGING_UVsur',                    &
     &   'Nudging surface currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NUDGING_UVsur,'
#endif
#ifdef NS_PERIODIC
      IF (Master) WRITE (stdout,20) 'NS_PERIODIC',                      &
     &   'North-South periodic boundaries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NS_PERIODIC,'
      FS_north=FS_north+1
      FS_south=FS_south+1
      M2_north=M2_north+1
      M2_south=M2_south+1
      M3_north=M3_north+1
      M3_south=M3_south+1
# ifdef WEC_VF
      M2s_north=M2s_north+1
      M2s_south=M2s_south+1
      M3s_north=M3s_north+1
      M3s_south=M3s_south+1
# endif
      K_north=K_north+1
      K_south=K_south+1
      T_north=T_north+1
      T_south=T_south+1
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
#ifdef OUT_DOUBLE
      IF (Master) WRITE (stdout,20) 'OUT_DOUBLE',                       &
     &   'Double precision output fields in NetCDF files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' OUT_DOUBLE,'
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
#if defined PARALLEL_IO && defined DISTRIBUTE
      IF (Master) WRITE (stdout,20) 'PARALLEL_IO',                      &
     &   'Parallel I/O processing with NetCDF-4/HDF5 libraries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' PERFECT_RESTART,'
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
#if defined Q_PSOURCE
      IF (Master) WRITE (stdout,20) 'Q_PSOURCE',                        &
     &   'Mass point sources, volume influx.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' Q_PSOURCE,'
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
#ifdef REFINED_GRID
      IF (Master) WRITE (stdout,20) 'REFINED_GRID',                     &
     &   'Refined grids option selected.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' REFINED_GRID,'
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
#if !(defined PJ_GRADPQ4 || defined PJ_GRADPQ2 || defined PJ_GRADP || \
      defined DJ_GRADPS) && defined RHO_SURF   && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'RHO_SURF',                         &
     &   'Include difference between rho0 and surface density.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' RHO_SURF,'
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
#ifdef SOUTH_FSCLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_FSCLAMPED',                  &
     &   'Southern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSCLAMPED,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSCHAPMAN
      IF (Master) WRITE (stdout,20) 'SOUTH_FSCHAPMAN',                  &
     &   'Southern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSCHAPMAN,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSGRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_FSGRADIENT',                 &
     &   'Southern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_FSGRADIENT,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSRADIATION
# ifdef SOUTH_FSNUDGING
      IF (Master) WRITE (stdout,20) 'SOUTH_FSNUDGING',                  &
     &   'Southern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'SOUTH_FSRADIATION',                &
     &   'Southern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_FSRADIATION,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_M2CLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_M2CLAMPED',                  &
     &   'Southern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2CLAMPED,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2FLATHER
      IF (Master) WRITE (stdout,20) 'SOUTH_M2FLATHER',                  &
     &   'Southern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2FLATHER,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2GRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_M2GRADIENT',                 &
     &   'Southern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M2GRADIENT,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2RADIATION
# ifdef SOUTH_M2NUDGING
      IF (Master) WRITE (stdout,20) 'SOUTH_M2NUDGING',                  &
     &   'Southern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'SOUTH_M2RADIATION',                &
     &   'Southern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M2RADIATION,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2REDUCED
      IF (Master) WRITE (stdout,20) 'SOUTH_M2REDUCED',                  &
     &   'Southern edge, 2D momentum, reduced-physics condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2REDUCED,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M3CLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_M3CLAMPED',                  &
     &   'Southern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M3CLAMPED,'
      M3_south=M3_south+1
#endif
#ifdef SOUTH_M3GRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_M3GRADIENT',                 &
     &   'Southern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M3GRADIENT,'
      M3_south=M3_south+1
#endif
#ifdef SOUTH_M3RADIATION
# ifdef SOUTH_M3NUDGING
      IF (Master) WRITE (stdout,20) 'SOUTH_M3NUDGING',                  &
     &   'Southern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M3NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'SOUTH_M3RADIATION',                &
     &   'Southern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M3RADIATION,'
      M3_south=M3_south+1
#endif
#ifdef SOUTH_M2SCLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_M2SCLAMPED',                 &
     &   'Southern edge, 2D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M2SCLAMPED,'
      M2s_south=M2s_south+1
#endif
#ifdef SOUTH_M2SGRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_M2SGRADIENT',                &
     &   'Southern edge, 2D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M2SGRADIENT,'
      M2s_south=M2s_south+1
#endif
#ifdef SOUTH_M2SRADIATION
      IF (Master) WRITE (stdout,20) 'SOUTH_M2SRADIATION',               &
     &   'Southern edge, 2D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' SOUTH_M2SRADIATION,'
      M2s_south=M2s_south+1
#endif
#ifdef SOUTH_M3SCLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_M3SCLAMPED',                 &
     &   'Southern edge, 3D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M3SCLAMPED,'
      M3s_south=M3s_south+1
#endif
#ifdef SOUTH_M3SGRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_M3SGRADIENT',                &
     &   'Southern edge, 3D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M3SGRADIENT,'
      M3s_south=M3s_south+1
#endif
#ifdef SOUTH_M3SRADIATION
      IF (Master) WRITE (stdout,20) 'SOUTH_M3SRADIATION',               &
     &   'Southern edge, 3D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' SOUTH_M3SRADIATION,'
      M3s_south=M3s_south+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef SOUTH_KGRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_KGRADIENT',                  &
     &   'Southern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_KGRADIENT,'
      K_south=K_south+1
# endif
# ifdef SOUTH_KRADIATION
      IF (Master) WRITE (stdout,20) 'SOUTH_KRADIATION',                 &
     &   'Southern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_KRADIATION,'
      K_south=K_south+1
# endif
#endif
#ifdef SOUTH_TCLAMPED
      IF (Master) WRITE (stdout,20) 'SOUTH_TCLAMPED',                   &
     &   'Southern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_TCLAMPED,'
      T_south=T_south+1
#endif
#ifdef SOUTH_TGRADIENT
      IF (Master) WRITE (stdout,20) 'SOUTH_TGRADIENT',                  &
     &   'Southern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_TGRADIENT,'
      T_south=T_south+1
#endif
#ifdef SOUTH_TRADIATION
# ifdef SOUTH_TNUDGING
      IF (Master) WRITE (stdout,20) 'SOUTH_TNUDGING',                   &
     &   'Southern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SOUTH_TNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'SOUTH_TRADIATION',                 &
     &   'Southern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_TRADIATION,'
      T_south=T_south+1
#endif
#ifdef SOUTH_VOLCONS
      IF (Master) WRITE (stdout,20) 'SOUTH_VOLCONS',                    &
     &   'Southern edge, enforce mass conservation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SOUTH_VOLCONS,'
#endif
#ifdef SOUTHERN_WALL
      IF (Master) WRITE (stdout,20) 'SOUTHERN_WALL',                    &
     &   'Wall boundary at Southern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SOUTHERN_WALL,'
#endif
#ifdef SO_SEMI
      IF (Master) WRITE (stdout,20) 'SO_SEMI',                          &
     &   'Stochastic Optimals, Semi-norm Estimation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SO_SEMI,'
      idriver=idriver+1
#endif
#ifdef SO_TRACE
      IF (Master) WRITE (stdout,20) 'SO_TRACE',                         &
     &   'Stochastic Optimals, Randomized Trace Estimation.'
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
#ifdef SPLINES
      IF (Master) WRITE (stdout,20) 'SPLINES',                          &
     &   'Conservative parabolic spline reconstruction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SPLINES,'
#endif
#if defined SPLINES_VCONV && defined FOUR_DVAR && defined VCONVOLUTION
      IF (Master) WRITE (stdout,20) 'IMPLICIT_VCONV',                   &
     &   'Implicit, Parabolic Splines Vertical Convolution Algorithm.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SPLINES_VCONV,'
#endif
#ifdef SPHERICAL
      IF (Master) WRITE (stdout,20) 'SPHERICAL',                        &
     &   'Spherical grid configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SPHERICAL,'
#endif
#ifdef SPONGE
      IF (Master) WRITE (stdout,20) 'SPONGE',                           &
     &   'Enhanced horizontal mixing in the sponge areas.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' SPONGE,'
#endif
#if defined SSH_TIDES
      IF (Master) WRITE (stdout,20) 'SSH_TIDES',                        &
     &   'Add tidal elevation to SSH climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SSH_TIDES,'
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
     &   'Nonconservative computation of alkalinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TALK_PROGNOSTIC,'
# else
      IF (Master) WRITE (stdout,20) '!TALK_PROGNOSTIC',                 &
     &   'Diagnostic computation of alkalinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !TALK_PROGNOSTIC,'
# endif
#endif
#ifdef TANGENT
      IF (Master) WRITE (stdout,20) 'TANGENT',                          &
     &   'Tangent Linear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TANGENT,'
#endif
#if defined TCLIMATOLOGY && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TCLIMATOLOGY',                     &
     &   'Processing tracer climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLIMATOLOGY,'
#endif
#if defined TCLM_NUDGING && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TCLM_NUDGING',                     &
     &   'Nudging toward tracer climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLM_NUDGING,'
#endif
#ifdef TLM_DRIVER
      IF (Master) WRITE (stdout,20) 'TLM_DRIVER',                       &
     &   'Generic tangent linear model driver.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' TLM_DRIVER,'
      idriver=idriver+1
#endif
#ifdef THREE_GHOST
      IF (Master) WRITE (stdout,20) 'THREE_GHOST',                      &
     &   'Using three Ghost Points in halo regions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' THREE_GHOST,'
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
# ifdef TS_MPDATA
      IF (Master) WRITE (stdout,20) 'TS_MPDATA',                        &
     &   'Recursive flux corrected MPDATA 3D advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' TS_MPDATA,'
      itrcHadv=itrcHadv+1
      itrcVadv=itrcVadv+1
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
       defined TS_U3HADVECTION)
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
       defined TS_SVADVECTION)
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
#if defined TS_PSOURCE && defined SOLVE3D
      IF (Master) WRITE (stdout,20) 'TS_PSOURCE',                       &
     &   'Tracers point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' TS_PSOURCE,'
#endif
#if defined TS_SMAGORINSKY && (defined TS_DIF2 || defined TS_DIF4)
      IF (Master) WRITE (stdout,20) 'TS_SMAGORINSKY',                   &
     &   'Smagorinksy-like time-dependent diffusion coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' TS_SMAGORINSKY,'
#endif
#ifdef UV_ADV
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
#endif
#if defined UV_BAROCLINIC && \
   (defined ASSIMILATION_UVsur || defined ASSIMILATION_UV)
      IF (Master) WRITE (stdout,20) 'UV_BAROCLINIC',                    &
     &   'Assimilation of only baroclinic currents.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_BAROCLINIC,'
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
#ifdef UV_PSOURCE
      IF (Master) WRITE (stdout,20) 'UV_PSOURCE',                       &
     &   'Mass point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' UV_PSOURCE,'
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
#if defined UV_SMAGORINSKY && (defined UV_DIF2 || defined UV_DIF4)
      IF (Master) WRITE (stdout,20) 'UV_SMAGORINSKY',                   &
     &   'Smagorinksy-like time-dependent viscosity coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_SMAGORINSKY,'
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
#ifdef WEST_FSCLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_FSCLAMPED',                   &
     &   'Western edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSCLAMPED,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSCHAPMAN
      IF (Master) WRITE (stdout,20) 'WEST_FSCHAPMAN',                   &
     &   'Western edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSCHAPMAN,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSGRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_FSGRADIENT',                  &
     &   'Western edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_FSGRADIENT,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSRADIATION
# ifdef WEST_FSNUDGING
      IF (Master) WRITE (stdout,20) 'WEST_FSNUDGING',                   &
     &   'Western edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'WEST_FSRADIATION',                 &
     &   'Western edge, free-sruface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_FSRADIATION,'
      FS_west=FS_west+1
#endif
#ifdef WEST_M2CLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_M2CLAMPED',                   &
     &   'Western edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2CLAMPED,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2FLATHER
      IF (Master) WRITE (stdout,20) 'WEST_M2FLATHER',                   &
     &   'Western edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2FLATHER,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2GRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_M2GRADIENT',                  &
     &   'Western edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M2GRADIENT,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2RADIATION
# ifdef WEST_M2NUDGING
      IF (Master) WRITE (stdout,20) 'WEST_M2NUDGING',                   &
     &   'Western edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'WEST_M2RADIATION',                 &
     &   'Western edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M2RADIATION,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2REDUCED
      IF (Master) WRITE (stdout,20) 'WEST_M2REDUCED',                   &
     &   'Western edge, 2D momentum, reduced-physics condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2REDUCED,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M3CLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_M3CLAMPED',                   &
     &   'Western edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M3CLAMPED,'
      M3_west=M3_west+1
#endif
#ifdef WEST_M3GRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_M3GRADIENT',                  &
     &   'Western edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M3GRADIENT,'
      M3_west=M3_west+1
#endif
#ifdef WEST_M3RADIATION
# ifdef WEST_M3NUDGING
      IF (Master) WRITE (stdout,20) 'WEST_M3NUDGING',                   &
     &   'Western edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M3NUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'WEST_M3RADIATION',                 &
     &   'Western edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M3RADIATION,'
      M3_west=M3_west+1
#endif
#ifdef WEST_M2SCLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_M2SCLAMPED',                  &
     &   'Westtern edge, 2D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M2SCLAMPED,'
      M2s_west=M2s_west+1
#endif
#ifdef WEST_M2SGRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_M2SGRADIENT',                 &
     &   'Western edge, 2D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M2SGRADIENT,'
      M2s_west=M2s_west+1
#endif
#ifdef WEST_M2SRADIATION
      IF (Master) WRITE (stdout,20) 'WEST_M2SRADIATION',                &
     &   'Western edge, 2D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' WEST_M2SRADIATION,'
      M2s_west=M2s_west+1
#endif
#ifdef WEST_M3SCLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_M3SCLAMPED',                  &
     &   'Western edge, 3D stokes, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M3SCLAMPED,'
      M3s_west=M3s_west+1
#endif
#ifdef WEST_M3SGRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_M3SGRADIENT',                 &
     &   'Western edge, 3D stokes, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M3SGRADIENT,'
      M3s_west=M3s_west+1
#endif
#ifdef WEST_M3SRADIATION
      IF (Master) WRITE (stdout,20) 'WEST_M3SRADIATION',                 &
     &   'Western edge, 3D stokes, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' WEST_M3SRADIATION,'
      M3s_west=M3s_west+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef WEST_KGRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_KGRADIENT',                   &
     &   'Western edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_KGRADIENT,'
      K_west=K_west+1
# endif
# ifdef WEST_KRADIATION
      IF (Master) WRITE (stdout,20) 'WEST_KRADIATION',                  &
     &   'Western edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_KRADIATION,'
      K_west=K_west+1
# endif
#endif
#ifdef WEST_TCLAMPED
      IF (Master) WRITE (stdout,20) 'WEST_TCLAMPED',                    &
     &   'Western edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_TCLAMPED,'
      T_west=T_west+1
#endif
#ifdef WEST_TGRADIENT
      IF (Master) WRITE (stdout,20) 'WEST_TGRADIENT',                   &
     &   'Western edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_TGRADIENT,'
      T_west=T_west+1
#endif
#ifdef WEST_TRADIATION
# ifdef WEST_TNUDGING
      IF (Master) WRITE (stdout,20) 'WEST_TNUDGING',                    &
     &   'Western edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' WEST_TNUDGING,'
# endif
      IF (Master) WRITE (stdout,20) 'WEST_TRADIATION',                  &
     &   'Western edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_TRADIATION,'
      T_west=T_west+1
#endif
#ifdef WEST_VOLCONS
      IF (Master) WRITE (stdout,20) 'WEST_VOLCONS',                     &
     &   'Western edge, enforce mass conservation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' WEST_VOLCONS,'
#endif
#ifdef WESTERN_WALL
      IF (Master) WRITE (stdout,20) 'WESTERN_WALL',                     &
     &   'Wall boundary at Western edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' WESTERN_WALL,'
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
#ifdef ZCLIMATOLOGY
      IF (Master) WRITE (stdout,20) 'ZCLIMATOLOGY',                     &
     &   'Processing sea surface height climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ZCLIMATOLOGY,'
#endif
#ifdef ZCLM_NUDGING
      IF (Master) WRITE (stdout,20) 'ZCLM_NUDGING',                     &
     &   'Nudging toward sea surface height climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ZCLM_NUDGING,'
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
      Coptions(is:is+14)=' ZOS_HSIG,'
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
!
!  Stop if more than one boundary conditions option for free-surface
!  is selected at each boundary edge.
!
      IF (Master.and.(FS_east.gt.1)) THEN
        WRITE (stdout,80) 'free-surface','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_west.gt.1)) THEN
        WRITE (stdout,80) 'free-surface','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_north.gt.1)) THEN
        WRITE (stdout,80) 'free-surface','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_south.gt.1)) THEN
        WRITE (stdout,80) 'free-surface','southern edge'
        exit_flag=5
      END IF
  80  FORMAT (/,' CHECKDEFS - only one boundary condition option for',  &
     &        1x,a,/,13x,'is allowed at the',1x,a,'.')
!
!  Stop if more than one boundary conditions option for 2D momentum
!  is selected at each boundary edge.
!
      IF (Master.and.(M2_east.gt.1)) THEN
        WRITE (stdout,80) '2D momentum','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_west.gt.1)) THEN
        WRITE (stdout,80) '2D momentum','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_north.gt.1)) THEN
        WRITE (stdout,80) '2D momentum','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_south.gt.1)) THEN
        WRITE (stdout,80) '2D momentum','southern edge'
        exit_flag=5
      END IF
#ifdef WEC_VF
!
!  Stop if more than one boundary conditions option for 2D stokes
!  is selected at each boundary edge.
!
      IF (Master.and.(M2s_east.gt.1)) THEN
        WRITE (stdout,80) '2D stokes','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2s_west.gt.1)) THEN
        WRITE (stdout,80) '2D stokes','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2s_north.gt.1)) THEN
        WRITE (stdout,80) '2D stokes','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2s_south.gt.1)) THEN
        WRITE (stdout,80) '2D stokes','southern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2s_east.eq.0)) THEN
        WRITE (stdout,20) 'EAST_M2SGRADIENT',                           &
     &  'Eastern edge, 2D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+28)=' EAST_M2SGRADIENT,'
      END IF
      IF (Master.and.(M2s_west.eq.0)) THEN
        WRITE (stdout,20) 'WEST_M2SGRADIENT',                           &
     &  'Western edge, 2D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+28)=' WEST_M2SGRADIENT,'
      END IF
      IF (Master.and.(M2s_north.eq.0)) THEN
        WRITE (stdout,20) 'NORTH_M2SGRADIENT',                          &
     &  'Northern edge, 2D stokes, gradientcondition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+29)=' NORTH_M2SGRADIENT,'
      END IF
      IF (Master.and.(M2s_south.eq.0)) THEN
        WRITE (stdout,20) 'SOUTH_M2SGRADIENT',                          &
     &  'Southern edge, 2D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+29)=' SOUTH_M2SGRADIENT,'
      END IF
#endif
#ifdef SOLVE3D
!
!  Stop if more than one boundary conditions option for 3D momentum
!  is selected at each boundary edge.
!
      IF (Master.and.(M3_east.gt.1)) THEN
        WRITE (stdout,80) '3D momentum','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_west.gt.1)) THEN
        WRITE (stdout,80) '3D momentum','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_north.gt.1)) THEN
        WRITE (stdout,80) '3D momentum','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_south.gt.1)) THEN
        WRITE (stdout,80) '3D momentum','southern edge'
        exit_flag=5
      END IF
#ifdef WEC_VF
!
!  Stop if more than one boundary conditions option for 3D stokes
!  is selected at each boundary edge.
!
      IF (Master.and.(M3s_east.gt.1)) THEN
        WRITE (stdout,80) '3D stokes','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3s_west.gt.1)) THEN
        WRITE (stdout,80) '3D stokes','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3s_north.gt.1)) THEN
        WRITE (stdout,80) '3D stokes','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3s_south.gt.1)) THEN
        WRITE (stdout,80) '3D stokes','southern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3s_east.eq.0)) THEN
        WRITE (stdout,20) 'EAST_M3SGRADIENT',                           &
     &  'Eastern edge, 3D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+28)=' EAST_M3SGRADIENT,'
      END IF
      IF (Master.and.(M3s_west.eq.0)) THEN
        WRITE (stdout,20) 'WEST_M3SGRADIENT',                           &
     &  'Western edge, 3D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+28)=' WEST_M3SGRADIENT,'
      END IF
      IF (Master.and.(M3s_north.eq.0)) THEN
        WRITE (stdout,20) 'NORTH_M3SGRADIENT',                          &
     &  'Northern edge, 3D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+29)=' NORTH_M3SGRADIENT,'
      END IF
      IF (Master.and.(M3s_south.eq.0)) THEN
        WRITE (stdout,20) 'SOUTH_M3SGRADIENT',                          &
     &  'Southern edge, 3D stokes, gradient condition.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+29)=' SOUTH_M3SGRADIENT,'
      END IF
#endif
!
!  Stop if more than one boundary conditions option for tracers is
!  selected at each boundary edge.
!
      IF (Master.and.(T_east.gt.1)) THEN
        WRITE (stdout,80) 'tracers','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_west.gt.1)) THEN
        WRITE (stdout,80) 'tracers','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_north.gt.1)) THEN
        WRITE (stdout,80) 'tracers','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_south.gt.1)) THEN
        WRITE (stdout,80) 'tracers','southern edge'
        exit_flag=5
      END IF
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Stop if more than one boundary conditions option for tracers is
!  selected at each boundary edge.
!
      IF (Master.and.(K_east.gt.1)) THEN
        WRITE (stdout,80) 'TKE fields','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_west.gt.1)) THEN
        WRITE (stdout,80) 'TKE fields','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_north.gt.1)) THEN
        WRITE (stdout,80) 'TKE fields','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_south.gt.1)) THEN
        WRITE (stdout,80) 'TKE fields','southern edge'
        exit_flag=5
      END IF
# endif
#endif
#ifndef SOLVE3D
!
!  Stop it explicit time-step splitting on shallow water set-up.
!
      DO ng=1,Ngrids
        IF (Master.and.(ndtfast(ng).gt.1)) THEN
          WRITE (stdout,90)
 90       FORMAT (/,' CHECKDEFS - explicit time-step splitting is ',    &
     &              ' inconsistent.',                                   &
     &            /,13x,'Change parameter NDTFAST to unity.')
          exit_flag=5
        END IF
      END DO
#endif
#if defined ADJOINT && defined _OPENMP
      IF (Master) THEN
        WRITE (stdout,100)
 100    FORMAT (/,' CHECKDEFS - cannot run the adjoint model in',       &
     &            ' shared-memory mode.')
        exit_flag=5
      END IF
#endif
#if defined NEMURO && defined HOLLING_GRAZING && defined IVLEV_EXPLICIT
!
!  Stop if using more than one grazing algorithm.
!
      IF (Master) THEN
        WRITE (stdout,110)
 110    FORMAT (/,' CHECKDEFS - More than one grazing algorithm',       &
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
        WRITE (stdout,120)
 120    FORMAT (/,' CHECKDEFS - More than one vertical convolution',    &
     &            ' algorithm selected.')
        exit_flag=5
      END IF
# endif
#endif
#if defined MODEL_COUPLING && (defined ESMF_LIB && defined MCT_LIB)
      IF (Master) THEN
        WRITE (stdout,130)
 130    FORMAT (/,' CHECKDEFS - More than one coupling library',        &
     &            ' selected.')
        exit_flag=5
      END IF
#endif
#if defined BULK_FLUXES && defined ATM2OCN_FLUXES
      IF (Master) THEN
        WRITE (stdout,130)
 130    FORMAT (/,' CHECKDEFS - More than one surface flux method ',    &
     &            ' selected.')
        exit_flag=5
      END IF
#endif
#if defined POSTERIOR_ERROR_F && defined POSTERIOR_ERROR_I && \
    defined WEAK_CONSTRAINT
      IF (Master) THEN
        WRITE (stdout,140)
 140    FORMAT (/,' CHECKDEFS - Both initial and final posterior',      &
     &            ' error covariance analysis selected, choose one.')
        exit_flag=5
      END IF
#endif
#ifdef UV_ADV
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,150) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
# ifdef SOLVE3D
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,150) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
# endif
#endif
#ifdef SOLVE3D
      IF (Master.and.(itrcHadv.gt.1)) THEN
        WRITE (stdout,150) 'horizontal','tracers','itrcHadv =',itrcHadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadv.gt.1)) THEN
        WRITE (stdout,150) 'vertical','tracers','itrcVadv =',itrcVadv
        exit_flag=5
      END IF
# if defined TANGENT || defined TL_IOMS || defined ADJOINT
      IF (Master.and.(itrcHadvtl.gt.1)) THEN
        WRITE (stdout,150) 'TL/AD horizontal','tracers','itrcHadvtl =', &
     &                     itrcHadvtl
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadvtl.gt.1)) THEN
        WRITE (stdout,150) 'TL/AD vertical','tracers','itrcVadvtl =',   &
     &                     itrcVadvtl
        exit_flag=5
      END IF
# endif
#endif
 150  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
#if defined SOLVE3D && defined NEMURO
# if !(defined IVLEV_EXPLICIT || defined HOLLING_GRAZING)
      IF (Master) THEN
        WRITE (stdout,160)
 160    FORMAT (/,' CHECKDEFS - Need to choose a zooplankton grazing',  &
     &            ' option;',/,13x,'use either IVLEV_EXPLICIT or',      &
     &            ' implicit HOLLING_GRAZING.',/,13x,'The default',     &
     &            ' implicit IVLEV algorithm does not work well yet.')
        exit_flag=5
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE checkdefs
