#include "cppdefs.h"
      MODULE bulk_flux_mod
#if defined BULK_FLUXES && !defined CCSM_FLUXES
!
!svn $Id: bulk_flux.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.     !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for       !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere       !
!      Response Experiment, JGR, 101, 3747-3764.                       !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.        !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer          !
!      effects on sea surface temperature, JGR, 101, 1295-1308.        !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
!  EMINUSP option for equivalent salt fluxes added by Paul Goodman     !
!  (10/2004).                                                          !
!                                                                      !
!  Modified by Kate Hedstrom for COARE version 3.0 (03/2005).          !
!  Modified by Jim Edson to correct specific hunidities.               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Fairall et al., 2003: J. Climate, 16, 571-591.                   !
!                                                                      !
!     Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea    !
!     surface roughness on the height and steepness of the waves.      !
!     J. Phys. Oceanogr., 31, 572-590.                                 !
!                                                                      !
!     Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002:!
!     New evidence for a relation between wind stress and wave age     !
!     from measurements during ASGAMAGE. Bound.-Layer Meteor., 103,    !
!     409-438.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: bulk_flux, bulk_psiu, bulk_psit
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bulk_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# ifdef ICE_MODEL
      USE mod_ice
# endif
# ifdef EVAP_SHALLOW
      USE mod_coupling
# endif
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL bulk_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef ICE_MODEL
     &                     liold(ng),                                   &
     &                     linew(ng),                                   &
# endif
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# endif
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % beta,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % t,                               &
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
# ifdef CLOUDS
     &                     FORCES(ng) % cloud,                          &
# endif
# if defined COARE_TAYLOR_YELLAND || defined DRENNAN
     &                     FORCES(ng) % Hwave,                          &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
     defined DRENNAN
     &                     FORCES(ng) % Pwave_top,                      &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
     &                     FORCES(ng) % Lwavep,                         &
# endif
# ifdef RUNOFF
     &                     FORCES(ng) % runoff,                         &
# endif
# ifdef EVAP_SHALLOW
     &                     GRID(ng) % h,                                &
     &                     COUPLING(ng) % Zt_avg1,                      &
# endif
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
     &                     FORCES(ng) % stflx,                          &
# ifdef CICE_MODEL
     &                     FORCES(ng) % LW_down,                        &
     &                     FORCES(ng) % SW_down,                        &
# endif
# if defined ALBEDO && defined SHORTWAVE
     &                     FORCES(ng) % albedo,                         &
#  ifdef ICE_MODEL
     &                     FORCES(ng) % albedo_ice,                     &
#  endif
# endif
# if defined ALBEDO_CLOUD && defined SHORTWAVE && !defined BENCHMARK
     &                     FORCES(ng) % cawdir,                         &
# endif
# ifdef ICE_MODEL
     &                     ICE(ng) % ai,                                &
     &                     ICE(ng) % hi,                                &
#  ifdef ICE_THERMO
     &                     FORCES(ng) % qai_n,                          &
     &                     FORCES(ng) % qi_o_n,                         &
     &                     FORCES(ng) % qswi_n,                         &
     &                     ICE(ng) % hsn,                               &
     &                     ICE(ng) % tis,                               &
     &                     ICE(ng) % coef_ice_heat,                     &
     &                     ICE(ng) % rhs_ice_heat,                      &
     &                     FORCES(ng) % snow_n,                         &
#   ifdef SNOWFALL
     &                     FORCES(ng) % snow,                           &
#   endif
#  endif
     &                     FORCES(ng) % sustr_aw,                       &
     &                     FORCES(ng) % svstr_aw,                       &
     &                     FORCES(ng) % qao_n,                          &
     &                     FORCES(ng) % tau_aix_n,                      &
     &                     FORCES(ng) % tau_aiy_n,                      &
# endif
# ifdef EMINUSP
     &                     FORCES(ng) % EminusP,                        &
     &                     FORCES(ng) % evap,                           &
# endif
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif
      RETURN
      END SUBROUTINE bulk_flux
!
!***********************************************************************
      SUBROUTINE bulk_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liold, linew,                          &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef WET_DRY
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
     &                           alpha, beta, rho, t,                   &
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
# ifdef CLOUDS
     &                           cloud,                                 &
# endif
# if defined COARE_TAYLOR_YELLAND || defined DRENNAN
     &                           Hwave,                                 &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
     defined DRENNAN
     &                           Pwave_top,                             &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
     &                           Lwavep,                                &
# endif
# ifdef RUNOFF
     &                           runoff,                                &
# endif
# ifdef EVAP_SHALLOW
     &                           h, Zt_avg1,                            &
# endif
     &                           rain, lhflx, lrflx, shflx,             &
     &                           srflx, stflx,                          &
# ifdef CICE_MODEL
     &                           LW_down, SW_down,                      &
# endif
# if defined ALBEDO && defined SHORTWAVE
     &                           albedo,                                &
#  ifdef ICE_MODEL
     &                           albedo_ice,                            &
#  endif
# endif
# if defined ALBEDO_CLOUD && defined SHORTWAVE && !defined BENCHMARK
     &                           cawdir,                                &
# endif
# ifdef ICE_MODEL
     &                           ai, hi,                                &
#  ifdef ICE_THERMO
     &                           qai_n, qi_o_n, qswi_n,                 &
     &                           hsn, tis, coef_ice_heat,               &
     &                           rhs_ice_heat, snow_n,                  &
#   ifdef SNOWFALL
     &                           snow,                                  &
#   endif
#  endif
     &                           sustr_aw, svstr_aw, qao_n,             &
     &                           tau_aix_n, tau_aiy_n,                  &
# endif
# ifdef EMINUSP
     &                           EminusP, evap,                         &
# endif
     &                           sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef ICE_MODEL
      USE mod_ice
# endif
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
# ifdef ICE_MODEL
      integer, intent(in) :: liold
      integer, intent(in) :: linew
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
# if defined COARE_TAYLOR_YELLAND || defined DRENNAN
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
      real(r8), intent(in) :: Lwavep(LBi:,LBj:)
# endif
      real(r8), intent(in) :: rain(LBi:,LBj:)
#  ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:,LBj:)
#  endif
#  ifdef EVAP_SHALLOW
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  endif

      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
# ifdef CICE_MODEL
      real(r8), intent(out) :: LW_down(LBi:,LBj:)
      real(r8), intent(out) :: SW_down(LBi:,LBj:)
# endif
# if defined ALBEDO && defined SHORTWAVE
      real(r8), intent(in) :: albedo(LBi:,LBj:)
#  ifdef ICE_MODEL
      real(r8), intent(in) :: albedo_ice(LBi:,LBj:)
#  endif
# endif
# if defined ALBEDO_CLOUD && defined SHORTWAVE && !defined BENCHMARK
      real(r8), intent(inout) :: cawdir(LBi:,LBj:)
#  endif

#  ifdef ICE_MODEL
      real(r8), intent(out) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(out) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(out) :: qao_n(LBi:,LBj:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qai_n(LBi:,LBj:)
      real(r8), intent(out) :: qi_o_n(LBi:,LBj:)
      real(r8), intent(out) :: qswi_n(LBi:,LBj:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
      real(r8), intent(out) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: rhs_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: snow_n(LBi:,LBj:)
#    ifdef SNOWFALL
      real(r8), intent(in) :: snow(LBi:,LBj:)
#    endif
#   endif
      real(r8), intent(out) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_aiy_n(LBi:,LBj:)
#   if defined ICE_BIO
      real(r8), intent(in) :: IcePhL(LBi:,LBj:,:)
      integer, intent(inout) :: IceLog(LBi:,LBj:,:)
#   endif
#  endif

#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:,LBj:)
      real(r8), intent(out) :: evap(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
# if defined COARE_TAYLOR_YELLAND || \
      defined DRENNAN
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
      real(r8), intent(in) :: Lwavep(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
#  ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:UBi,LBj:UBj)
#  endif
#  ifdef EVAP_SHALLOW
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# ifdef CICE_MODEL
      real(r8), intent(out) :: LW_down(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SW_down(LBi:UBi,LBj:UBj)
# endif
# if defined ALBEDO && defined SHORTWAVE
      real(r8), intent(inout) :: albedo(LBi:UBi,LBj:UBj)
#  ifdef ICE_MODEL
      real(r8), intent(inout) :: albedo_ice(LBi:UBi,LBj:UBj)
#  endif
# endif
# if defined ALBEDO_CLOUD && defined SHORTWAVE && !defined BENCHMARK
      real(r8), intent(inout) :: cawdir(LBi:UBi,LBj:UBj)
#  endif

#  ifdef ICE_MODEL
      real(r8), intent(out) :: sustr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qao_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qai_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qi_o_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qswi_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: coef_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhs_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: snow_n(LBi:UBi,LBj:UBj)
#    ifdef SNOWFALL
      real(r8), intent(in) :: snow(LBi:UBi,LBj:UBj)
#    endif
#   endif
      real(r8), intent(out) :: tau_aix_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tau_aiy_n(LBi:UBi,LBj:UBj)
#   if defined ICE_BIO
      real(r8), intent(in) :: IcePhL(LBi:UBi,LBj:UBj,2)
      integer, intent(inout) :: IceLog(LBi:UBi,LBj:UBj,2)
#   endif
#  endif

#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k, listp
      integer :: iday, month, year
      real(r8) :: hour, yday
      integer, parameter :: IterMax = 3

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8), parameter :: Ch_i = 1.75e-3
      real(r8), parameter :: Ce_i = 1.75e-3
! ratio of molecular weight of water to dry air
      real(r8), parameter :: epsilon = 0.622_r8

      real(r8) :: Bf, Cd, Ch, Hl, Hlw, Hscale, Hscale2, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL
# ifdef ICE_BULK_FLUXES
      real(r8) :: qswi, qlwi, qlh_i, qsh_i
      real(r8) :: le_i, dq_i,fqlat1, sfc_temp, slp, Qsati
      real(r8) :: vap_p_i
      real(r8), parameter :: Io_frac = 0.17_r8   ! from MU
# endif

      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(r8) :: twopi_inv, wet_bulb
# if defined SHORTWAVE && defined ALBEDO_CLOUD
      real(r8), parameter :: albw_d=0.065_r8
      real(r8) :: cawdif
# endif
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
# endif
#ifdef ICE_BOX
      real(r8), parameter :: sensible_MU(12) =                          &
     &        (/ 19.1, 12.3, 11.6, 4.7, -7.3, -6.3,                     &
     &           -4.8, -6.5, -2.7, 1.6, 9.0, 12.8     /)
      real(r8), parameter :: latent_MU(12) =                            &
     &        (/ 0.0, -0.32, -0.48, -1.45, -7.4, -11.3,                &
     &           -10.3, -10.7, -6.3, -3.1, -.16, -.16   /)
#endif

      real(r8), dimension(IminS:ImaxS) :: CC
      real(r8), dimension(IminS:ImaxS) :: Cd10
      real(r8), dimension(IminS:ImaxS) :: Ch10
      real(r8), dimension(IminS:ImaxS) :: Ct10
      real(r8), dimension(IminS:ImaxS) :: charn
      real(r8), dimension(IminS:ImaxS) :: Ct
      real(r8), dimension(IminS:ImaxS) :: Cwave
      real(r8), dimension(IminS:ImaxS) :: Cwet
      real(r8), dimension(IminS:ImaxS) :: delQ
      real(r8), dimension(IminS:ImaxS) :: delQc
      real(r8), dimension(IminS:ImaxS) :: delT
      real(r8), dimension(IminS:ImaxS) :: delTc
      real(r8), dimension(IminS:ImaxS) :: delW
      real(r8), dimension(IminS:ImaxS) :: L
      real(r8), dimension(IminS:ImaxS) :: L10
      real(r8), dimension(IminS:ImaxS) :: Q
      real(r8), dimension(IminS:ImaxS) :: Qair
      real(r8), dimension(IminS:ImaxS) :: Qpsi
      real(r8), dimension(IminS:ImaxS) :: Qsea
      real(r8), dimension(IminS:ImaxS) :: Qstar
      real(r8), dimension(IminS:ImaxS) :: rhoAir
      real(r8), dimension(IminS:ImaxS) :: rhoSea
      real(r8), dimension(IminS:ImaxS) :: Ri
      real(r8), dimension(IminS:ImaxS) :: Ribcu
      real(r8), dimension(IminS:ImaxS) :: Rr
      real(r8), dimension(IminS:ImaxS) :: Scff
      real(r8), dimension(IminS:ImaxS) :: TairC
      real(r8), dimension(IminS:ImaxS) :: TairK
      real(r8), dimension(IminS:ImaxS) :: Tcff
      real(r8), dimension(IminS:ImaxS) :: Tpsi
      real(r8), dimension(IminS:ImaxS) :: TseaC
      real(r8), dimension(IminS:ImaxS) :: TseaK
      real(r8), dimension(IminS:ImaxS) :: Tstar
      real(r8), dimension(IminS:ImaxS) :: u10
      real(r8), dimension(IminS:ImaxS) :: VisAir
      real(r8), dimension(IminS:ImaxS) :: WaveLength
      real(r8), dimension(IminS:ImaxS) :: Wgus
      real(r8), dimension(IminS:ImaxS) :: Wmag
      real(r8), dimension(IminS:ImaxS) :: Wpsi
      real(r8), dimension(IminS:ImaxS) :: Wstar
      real(r8), dimension(IminS:ImaxS) :: Zetu
      real(r8), dimension(IminS:ImaxS) :: Zo10
      real(r8), dimension(IminS:ImaxS) :: ZoT10
      real(r8), dimension(IminS:ImaxS) :: ZoL
      real(r8), dimension(IminS:ImaxS) :: ZoQ
      real(r8), dimension(IminS:ImaxS) :: ZoT
      real(r8), dimension(IminS:ImaxS) :: ZoW
      real(r8), dimension(IminS:ImaxS) :: ZWoL

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy
#ifdef ICE_THERMO
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Coef_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: RHS_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qai
#endif
#ifdef ICE_MODEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
#endif
#ifdef EVAP_SHALLOW
      real(r8), parameter :: coef_shallow = 5.28_r8
!      real(r8), parameter :: coef_shallow = 4.69_r8
      real(r8) :: cffe
#endif

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
      twopi_inv=0.5_r8/pi
#if defined ICE_MODEL
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        listp = liold
      ELSE
        listp = linew
      END IF
#endif

#ifdef ICE_BOX
      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
#endif

      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))+eps
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
# if defined ICE_MODEL && defined ICE_THERMO
          qswi_n(i,j) = SRad(i,j)
# endif
# if defined ALBEDO && defined SHORTWAVE
#   if defined ALBEDO_CLOUD && !defined BENCHMARK
          SRad(i,j) = SRad(i,j)*cawdir(i,j)
#  else
          SRad(i,j) = SRad(i,j)*(1.0_r8-albedo(i,j))
#  endif
# endif
          Tcff(i)=alpha(i,j)
          Scff(i)=beta(i,j)
!
!  Initialize.
!
          delTc(i)=0.0_r8
          delQc(i)=0.0_r8
          LHeat(i,j)=lhflx(i,j)*Hscale
          SHeat(i,j)=shflx(i,j)*Hscale
          Taur=0.0_r8
          Taux(i,j)=0.0_r8
          Tauy(i,j)=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------

# if defined LONGWAVE
!
!  Use Berliand (1952) formula to calculate net longwave radiation.
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
!  is assumed constant, but it is a function of latitude varying from
!  1.0 at poles to 0.5 at the Equator).
!
#  ifdef SPECIFIC_HUMIDITY
!  specific humidity in units of kg/kg
          vap_p=PairM*RH/(1.0_r8+0.378_r8*RH)             !WPB
#  else
      IF(RH.lt.2.0_r8) THEN                               !WPB
          cff=(0.7859_r8+0.03477_r8*TairC(i))/                          &
     &        (1.0_r8+0.00412_r8*TairC(i))
          e_sat=10.0_r8**cff   ! saturation vapor pressure (hPa or mbar)
          vap_p=e_sat*RH       ! water vapor pressure (hPa or mbar)
      ELSE                                                !WPB
          vap_p=0.001_r8*PairM*RH/(1.0_r8+0.000378_r8*RH) !WPB
      ENDIF                                               !WPB
#  endif
          cff2=TairK(i)*TairK(i)*TairK(i)
          cff1=cff2*TairK(i)
          LRad(i,j)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_r8-0.05_r8*SQRT(vap_p))*                &
!    &              (cff1*(0.39_r8-0.05_r8*SQRT(0.01_r8*vap_p))*        &
     &                    (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))+     &
     &               cff2*4.0_r8*(TseaK(i)-TairK(i)))

# elif defined LONGWAVE_OUT
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          LRad(i,j)=lrflx(i,j)*Hscale-                                  &
     &              emmiss*StefBo*TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)

# else
          LRad(i,j)=lrflx(i,j)*Hscale
# endif
# ifdef MASKING
          LRad(i,j)=LRad(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          LRad(i,j)=LRad(i,j)*rmask_wet(i,j)
# endif
!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair(i) is the saturation specific humidity at Tair
!                 Q(i) is the actual specific humidity
!              Qsea(i) is the saturation specific humidity at Tsea
!
!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007_r8+3.46E-6_r8*PairM(mb))*6.1121_r8*
!                  EXP(17.502_r8*TairC(C)/(240.97_r8+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197_r8*(E(mb)/(PairM(mb)-0.378_r8*E(mb)))
!
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TairC(i)/(240.97_r8+TairC(i)))
!
!  Compute specific humidity at Saturation, Qair (kg/kg).
!
          Qair(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff+eps))
!
!  Compute specific humidity, Q (kg/kg).
!
#if defined CORE_FORCING || defined SPECIFIC_HUMIDITY
! Incoming humidity is specific humidity in (kg/kg)
          Q(i)=RH
#else
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff+eps)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
#endif
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
          cff=cff*0.98_r8
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
          Qsea(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir(i)=PairM*100.0_r8/(blk_Rgas*TairK(i)*                  &
     &                              (1.0_r8+0.61_r8*Q(i)))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir(i)=1.326E-5_r8*                                        &
     &              (1.0_r8+TairC(i)*(6.542E-3_r8+TairC(i)*             &
     &               (8.301E-6_r8-4.84E-9_r8*TairC(i))))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv(i,j)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus(i)=0.5_r8
          delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          delQ(i)=Qsea(i)-Q(i)
          delT(i)=TseaC(i)-TairC(i)
!
!  Neutral coefficients.
!
          ZoW(i)=0.0001_r8
          u10(i)=delW(i)*LOG(10.0_r8/ZoW(i))/LOG(blk_ZW(ng)/ZoW(i))
          Wstar(i)=0.035_r8*u10(i)
          Zo10(i)=0.011_r8*Wstar(i)*Wstar(i)/g+                         &
     &            0.11_r8*VisAir(i)/Wstar(i)
          Cd10(i)=(vonKar/LOG(10.0_r8/Zo10(i)))**2
          Ch10(i)=0.00115_r8
          Ct10(i)=Ch10(i)/sqrt(Cd10(i))
          ZoT10(i)=10.0_r8/EXP(vonKar/Ct10(i))
          Cd=(vonKar/LOG(blk_ZW(ng)/Zo10(i)))**2
!
!  Compute Richardson number.
!
          Ct(i)=vonKar/LOG(blk_ZT(ng)/ZoT10(i))  ! T transfer coefficient
          CC(i)=vonKar*Ct(i)/Cd
          delTc(i)=0.0_r8
# ifdef COOL_SKIN
          delTc(i)=blk_dter
# endif
          Ribcu(i)=-blk_ZW(ng)/(blk_Zabl*0.004_r8*blk_beta**3)
          Ri(i)=-g*blk_ZW(ng)*((delT(i)-delTc(i))+                      &
     &                          0.61_r8*TairK(i)*delQ(i))/              &
     &          (TairK(i)*delW(i)*delW(i)+eps)
          IF (Ri(i).lt.0.0_r8) THEN
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+Ri(i)/Ribcu(i))       ! Unstable
          ELSE
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+3.0_r8*Ri(i)/CC(i))   ! Stable
          END IF
          L10(i)=blk_ZW(ng)/Zetu(i)
!
!  First guesses for Monon-Obukhov similarity scales.
!
          Wstar(i)=delW(i)*vonKar/(LOG(blk_ZW(ng)/Zo10(i))-             &
     &                             bulk_psiu(blk_ZW(ng)/L10(i),pi))
          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                          &
     &             (LOG(blk_ZT(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZT(ng)/L10(i),pi))
          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                          &
     &             (LOG(blk_ZQ(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZQ(ng)/L10(i),pi))
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          IF (delW(i).gt.18.0_r8) THEN
            charn(i)=0.018_r8
          ELSE IF ((10.0_r8.lt.delW(i)).and.(delW(i).le.18.0_r8)) THEN
            charn(i)=0.011_r8+                                          &
     &               0.125_r8*(0.018_r8-0.011_r8)*(delW(i)-10._r8)
          ELSE
            charn(i)=0.011_r8
          END IF
# if defined COARE_OOST || defined COARE_TAYLOR_YELLAND ||   \
     defined DRENNAN
#  if defined DEEPWATER_WAVES
          Cwave(i)=g*MAX(Pwave_top(i,j),eps)*twopi_inv
          WaveLength(i)=Cwave(i)*MAX(Pwave_top(i,j),eps)
#  else
          Cwave(i)=Lwavep(i,j)/MAX(Pwave_top(i,j),eps)
          WaveLength(i)=Lwavep(i,j)
#  endif
# endif
        END DO
!
!  Iterate until convergence. It usually converges within 3 iterations.
# if defined COARE_OOST || defined COARE_TAYLOR_YELLAND
!  Use wave info if we have it, two different options.
# endif
!
        DO Iter=1,IterMax
          DO i=Istr-1,IendR
# ifdef COARE_OOST
            ZoW(i)=(25.0_r8/pi)*WaveLength(i)*                          &
     &             (Wstar(i)/Cwave(i))**4.5_r8+                         &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# elif defined COARE_TAYLOR_YELLAND
            ZoW(i)=1200.0_r8*Hwave(i,j)*                                &
     &             (Hwave(i,j)/WaveLength(i))**4.5_r8+                  &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# elif defined DRENNAN
            ZoW(i)=(3.35_r8)*Hwave(i,j)*                                &
     &             MIN((Wstar(i)/Cwave(i)),0.1_r8)**3.4_r8+             &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# else
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6_r8)
            ZoT(i)=ZoQ(i)
            ZoL(i)=vonKar*g*blk_ZW(ng)*                                 &
     &             (Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                     &
     &                        0.61_r8*TairK(i)*Qstar(i))/               &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
            L(i)=blk_ZW(ng)/(ZoL(i)+eps)
!
!  Evaluate stability functions at Z/L.
!
            Wpsi(i)=bulk_psiu(ZoL(i),pi)
            Tpsi(i)=bulk_psit(blk_ZT(ng)/L(i),pi)
            Qpsi(i)=bulk_psit(blk_ZQ(ng)/L(i),pi)
# ifdef COOL_SKIN
            Cwet(i)=0.622_r8*Hlv(i,j)*Qsea(i)/                          &
     &              (blk_Rgas*TseaK(i)*TseaK(i))
            delQc(i)=Cwet(i)*delTc(i)
# endif
!
!  Compute wind scaling parameters, Wstar.
!
            Wstar(i)=MAX(eps,delW(i)*vonKar/                            &
     &               (LOG(blk_ZW(ng)/ZoW(i))-Wpsi(i)))
            Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
     &               (LOG(blk_ZT(ng)/ZoT(i))-Tpsi(i))
            Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
     &               (LOG(blk_ZQ(ng)/ZoQ(i))-Qpsi(i))
!
!  Compute gustiness in wind speed.
!
            Bf=-g/TairK(i)*                                             &
     &         Wstar(i)*(Tstar(i)+0.61_r8*TairK(i)*Qstar(i))
            IF (Bf.gt.0.0_r8) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
            ELSE
              Wgus(i)=0.2_r8
            END IF
            delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
# ifdef COOL_SKIN
!
!-----------------------------------------------------------------------
!  Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
            Clam=16.0_r8*g*blk_Cpw*(rhoSea(i)*blk_visw)**3/             &
     &           (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
            Hcool=0.001_r8
!
!  Backgound sensible and latent heat.
!
            Hsb=-rhoAir(i)*blk_Cpa*Wstar(i)*Tstar(i)
            Hlb=-rhoAir(i)*Hlv(i,j)*Wstar(i)*Qstar(i)
!
!  Mean absoption in cool-skin layer.
!
            Fc=0.065_r8+11.0_r8*Hcool-                                  &
     &         (1.0_r8-EXP(-Hcool*1250.0_r8))*6.6E-5_r8/Hcool
!
!  Total cooling at the interface.
!
            Qcool=LRad(i,j)+Hsb+Hlb-SRad(i,j)*Fc
            Qbouy=Tcff(i)*Qcool+Scff(i)*Hlb*blk_Cpw/Hlv(i,j)
!
!  Compute temperature and moisture change.
!
            IF ((Qcool.gt.0.0_r8).and.(Qbouy.gt.0.0_r8)) THEN
              lambd=6.0_r8/                                             &
     &              (1.0_r8+                                            &
     &               (Clam*Qbouy/(Wstar(i)+eps)**4)**0.75_r8)**r3
              Hcool=lambd*blk_visw/(SQRT(rhoAir(i)/rhoSea(i))*          &
     &                              Wstar(i)+eps)
              delTc(i)=Qcool*Hcool/blk_tcw
            ELSE
              delTc(i)=0.0_r8
            END IF
            delQc(i)=Cwet(i)*delTc(i)
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
        DO i=Istr-1,IendR
!
!  Compute transfer coefficients for momentum (Cd).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
          Ch=Wstar(i)*Tstar(i)/(-Wspeed*delT(i)+0.0098_r8*blk_ZT(ng))
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
          Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
          diffw=2.11E-5_r8*(TairK(i)/273.16_r8)**1.94_r8
          diffh=0.02411_r8*(1.0_r8+TairC(i)*                            &
     &                      (3.309E-3_r8-1.44E-6_r8*TairC(i)))/         &
     &          (rhoAir(i)*blk_Cpa)
          cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
          wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/        &
     &                                     (blk_Cpa*diffh))
          Hsr=rain(i,j)*wet_bulb*blk_Cpw*                               &
     &        ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
          SHeat(i,j)=(Hs+Hsr)
# ifdef WET_DRY
          SHeat(i,j)=SHeat(i,j)*rmask_wet(i,j)
# elif defined MASKING
          SHeat(i,j)=SHeat(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          SHeat(i,j)=SHeat(i,j)*rmask_wet(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
          Hl=-Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61_r8*Wstar(i)*Qstar(i)-                             &
     &          (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
          Hlw=rhoAir(i)*Hlv(i,j)*upvel*Q(i)
          LHeat(i,j)=(Hl+Hlw)
# ifdef ENHANCED_LHF_TRANSFER
          LHeat(i,j)=LHeat(i,j)*1.07_r8
# endif
# ifdef WET_DRY
          LHeat(i,j)=LHeat(i,j)*rmask_wet(i,j)
# elif defined MASKING
          LHeat(i,j)=LHeat(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          LHeat(i,j)=LHeat(i,j)*rmask_wet(i,j)
# endif
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
          Taur=0.85_r8*rain(i,j)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
          cff=rhoAir(i)*Cd*Wspeed
          Taux(i,j)=(cff*Uwind(i,j)+Taur*SIGN(1.0_r8,Uwind(i,j)))
# ifdef WET_DRY
          Taux(i,j)=Taux(i,j)*rmask_wet(i,j)
# elif defined MASKING
          Taux(i,j)=Taux(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          Taux(i,j)=Taux(i,j)*rmask_wet(i,j)
# endif
          Tauy(i,j)=(cff*Vwind(i,j)+Taur*SIGN(1.0_r8,Vwind(i,j)))
# ifdef WET_DRY
          Tauy(i,j)=Tauy(i,j)*rmask_wet(i,j)
# elif defined MASKING
          Tauy(i,j)=Tauy(i,j)*rmask(i,j)
# endif
# ifdef ICE_MODEL
          Taux_Ice(i,j) = rhoAir(i)*cdai(ng)*Uwind(i,j)*Wspeed
          Tauy_Ice(i,j) = rhoAir(i)*cdai(ng)*Vwind(i,j)*Wspeed

#  ifdef ICE_THERMO
!  Correct ocean net short-wave radiation for ice cover and convert to
!  kinematic units

          Coef_Ice(i,j) = 0._r8
          RHS_Ice(i,j) = 0._r8

! Sensible heat flux over ice

          IF (ai(i,j,listp).gt.min_a(ng)) THEN
            sfc_temp = tis(i,j) + 273.16_r8
          ELSE
            sfc_temp = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          ENDIF
#   ifdef ICE_BOX
          qsh_i = sensible_MU(month)
          RHS_Ice(i,j) = RHS_Ice(i,j) + qsh_i
#   else
          qsh_i = rhoAir(i)*blk_Cpa*Ch_i*delW(i)*                         &
     &            (TairK(i)+0.0098*blk_ZT(ng)-sfc_temp)
          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &          rhoAir(i)*blk_Cpa*Ch_i*delW(i)
          RHS_Ice(i,j) = RHS_Ice(i,j) +                                 &
     &        rhoAir(i)*blk_Cpa*Ch_i*delW(i)*(TairK(i)+0.0098*blk_ZT(ng))
#   endif

! Latent heat flux over ice
          le_i=2.834E+6_r8    ! Latent heat of sublimation

#   ifdef ICE_BOX
          qlh_i = latent_MU(month)
#   else
! Qsati is saturation specific humidity over the ice from
! Parkinson and Washington (1979)
          slp = Pair(i,j)*100._r8 !Convert back to Pascal from millibars
          cff = 611._r8*                                                &
     &      10._r8**(9.5_r8*(sfc_temp-273.16_r8)/(sfc_temp-7.66_r8))
!         Qsati=epsilon*cff/(slp-(1._r8-.62197_r8)*cff)
          Qsati=0.62197_r8*cff/(slp-(1._r8-.62197_r8)*cff)
          dq_i=Q(i)-Qsati
          fqlat1=rhoAir(i)*Ce_i*delW(i)
          qlh_i = fqlat1*le_i*dq_i
#   endif
          RHS_Ice(i,j) = RHS_Ice(i,j) + qlh_i
!
! Short-wave radiation to ice
#   if defined ALBEDO && defined SHORTWAVE
          qswi = (1._r8-albedo_ice(i,j))*qswi_n(i,j)
#   else
          qswi = qswi_n(i,j)
#   endif

! Already in W/m**2
          RHS_Ice(i,j) = RHS_Ice(i,j) + qswi


! Net long-wave radiation from ice
          cff = Q(i)/(1.0_r8+Q(i))
          vap_p_i = slp*cff/(0.62197 + cff)
#   ifdef LONGWAVE
          qlwi = emmiss*StefBo*(TairK(i)**4)*                           &
     &         (0.39_r8-0.005_r8*SQRT(vap_p_i))*                        &
     &                (1.0_r8-0.65_r8*cloud(i,j)) -                     &
     &          4._r8*StefBo*emmiss*TairK(i)**3*(sfc_temp-TairK(i))

          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &           4._r8*StefBo*emmiss*TairK(i)**3
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &         emmiss*StefBo*(TairK(i)**4)*(                            &
     &    (0.39_r8-0.005_r8*SQRT(vap_p_i))*(1.0_r8-0.65_r8*cloud(i,j))  &
     &    - 4.0_r8)
#   else
          qlwi = lrflx(i,j)*Hscale - emmiss*StefBo*sfc_temp**4
          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &            4._r8*StefBo*emmiss*sfc_temp**3
          RHS_Ice(i,j) = RHS_Ice(i,j) + lrflx(i,j)*Hscale +             &
     &            3._r8*StefBo*emmiss*sfc_temp**4
#   endif
! Correct for Kelvin to Celsius in RHS
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &            Coef_Ice(i,j)*273.16_r8

! Net heat flux from ice to atmosphere
          Qai(i,j) = -qsh_i - qlh_i - qswi - qlwi
#  endif
# endif
          !--- heat flux ---
# ifdef ICE_BOX
          SHeat(i,j) = sensible_MU(month)
          LHeat(i,j) = latent_MU(month)
# else
#  ifdef ICE_MODEL
          SHeat(i,j) = SHeat(i,j)*(1-ai(i,j,listp)) +                   &
     &                    ai(i,j,listp)*qsh_i
          LHeat(i,j) = LHeat(i,j)*(1-ai(i,j,listp)) +                   &
     &                    ai(i,j,listp)*qlh_i
#  endif
# endif
        END DO
      END DO
!
!=======================================================================
!  Compute surface net heat flux and surface wind stress.
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  calculated from the bulk formulations above are positive out of the
!  ocean.
!
!  For EMINUSP option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
!                       PREC = rain = kg/m2/s
!
!  To convert these rates to m/s divide by freshwater density, rhow.
!
!  Note that when the air is undersaturated in water vapor (Q < Qsea)
!  the model will evaporate and LHeat > 0:
!
!                   LHeat positive out of the ocean
!                    evap positive out of the ocean
!
!  Note that if evaporating, the salt flux is positive
!        and if     raining, the salt flux is negative
!
!  Note that fresh water flux is positive out of the ocean and the
!  salt flux (stflx(isalt)) is positive into the ocean. It is converted
!  to (psu m/s) for stflx(isalt) in "set_vbc.F" or "ice_mk.h". The E-P
!  value is saved in variable EminusP for I/O purposes.
!
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
# ifdef CICE_MODEL
          LW_down(i,j) = lrflx(i,j)*Hscale
          SW_down(i,j) = srflx(i,j)*Hscale
# endif
# if defined ALBEDO && defined SHORTWAVE
#   if defined ALBEDO_CLOUD && !defined BENCHMARK
          cawdif = 1._r8-albw_d
          srflx(i,j) = srflx(i,j)*(cawdir(i,j)*                         &
     &                 (1._r8-cloud(i,j))+cawdif*cloud(i,j))
#  else
          srflx(i,j) = srflx(i,j)*(1._r8-albedo(i,j))
#  endif
# endif
          lrflx(i,j) = LRad(i,j)*Hscale2
          lhflx(i,j) = -LHeat(i,j)*Hscale2
          shflx(i,j) = -SHeat(i,j)*Hscale2
          stflx(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                      &
     &                      lhflx(i,j)+shflx(i,j))
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
# endif
# ifdef EMINUSP
          evap(i,j)=LHeat(i,j)/(Hlv(i,j)+eps)
          stflx(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
#  ifdef ICE_MODEL
          evap(i,j) = (1.0_r8-ai(i,j,listp))*evap(i,j)
          stflx(i,j,isalt) = (1.0_r8-ai(i,j,listp))*stflx(i,j,isalt)
#  endif
#  ifdef EVAP_SHALLOW
          cffe = 1._r8/(h(i,j)+Zt_avg1(i,j))
          cffe = min(cffe,0.2_r8)
          IF (evap(i,j).GT.0._r8) THEN
             evap(i,j) = (1.0_r8+coef_shallow*cffe)*evap(i,j)
          END IF
#  endif
# endif
# ifdef WET_DRY
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask_wet(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#   endif
# elif defined MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#   endif
# endif
# ifdef RUNOFF
          stflx(i,j,isalt) = stflx(i,j,isalt) -                         &
     &                        cff*runoff(i,j)
#  ifdef WET_DRY
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#  elif defined MASKING
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#  endif
# endif

! Limit shortwave radiation to be used in computing penetrative
! radiation by presence of ice and biology

# if defined ICE_MODEL
! Limit shortwave radiation to be used in computing penetrative
! radiation by presence of ice alone
          ice_thick(i,j) = hi(i,j,listp)/(ai(i,j,listp)+eps)
          snow_thick(i,j) = hsn(i,j,listp)/(ai(i,j,listp)+eps)

          IF (ice_thick(i,j).le.0.1_r8) THEN
            cff1=ice_thick(i,j)*5._r8
          ELSE
            cff1=0.1_r8*5._r8+(ice_thick(i,j)-0.1_r8)*1._r8
          ENDIF

          cff1= cff1+snow_thick(i,j)*20._r8

#  ifdef ICE_BIO
          cff2=0.5*0.02_r8*0.005*IcePhL(i,j,listp)
          srflx(i,j)=(1.0_r8-ai(i,j,listp))*srflx(i,j)                  &
     &            +ai(i,j,listp)*srflx(i,j)*exp(-cff1-cff2)
#  elif !defined ICE_BOX
          srflx(i,j)=(1.0_r8-ai(i,j,listp))*srflx(i,j)                  &
     &            +ai(i,j,listp)*srflx(i,j)*exp(-cff1)
#  endif
          qao_n(i,j) = -stflx(i,j,itemp)*Hscale

! Net precipitation - evaporation in m/s
          IF (Tair(i,j).ge.0.0) THEN
            snow_n(i,j) = 0.0_r8
          ELSE
#  ifdef SNOWFALL
            snow_n(i,j) = snow(i,j)/rhosnow_dry(ng)
#  else
            snow_n(i,j) = rain(i,j)/rhosnow_dry(ng)
#  endif
          ENDIF
          coef_ice_heat(i,j) = Coef_Ice(i,j)
          rhs_ice_heat(i,j) = RHS_Ice(i,j)
          qai_n(i,j) = Qai(i,j)

# endif

# ifdef EMINUSP
          EminusP(i,j)=stflx(i,j,isalt)
# endif
        END DO
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      cff=0.5_r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# elif defined MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef ICE_MODEL
          sustr_aw(i,j) = sustr(i,j)
          tau_aix_n(i,j) = Taux_Ice(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# elif defined MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef ICE_MODEL
          svstr_aw(i,j) = svstr(i,j)
          tau_aiy_n(i,j) = Tauy_Ice(i,j)
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lrflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))
# ifdef CICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          LW_down)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          SW_down)
# endif
# ifdef ICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qao_n)
#  ifdef ICE_THERMO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qai_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qi_o_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          coef_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhs_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          snow_n)
#  endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_aix_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_aiy_n)
# endif
# if defined EMINUSP || (defined ICE_MODEL && defined ICE_THERMO)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          EminusP)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
# ifdef ICE_MODEL
        CALL exchange_u2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr_aw)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr_aw)
# endif
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflx(:,:,itemp))
#  ifdef CICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    LW_down, SW_down)
#  endif
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, qao_n, tau_aix_n, tau_aiy_n )
#   ifdef ICE_THERMO
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    qai_n, coef_ice_heat, rhs_ice_heat, snow_n )
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    qi_o_n)
#   endif
#  endif
#  if defined EMINUSP || (defined ICE_MODEL && defined ICE_THERMO)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap, EminusP,                                &
     &                    stflx(:,:,isalt))
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr_aw, svstr_aw)
#  endif
# endif

      RETURN
      END SUBROUTINE bulk_flux_tile

      FUNCTION bulk_psiu (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
!  Function result
!
      real(r8) :: bulk_psiu
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.25_r8
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))+                             &
     &       LOG(0.5_r8*(1.0_r8+x*x))-                                  &
     &       2.0_r8*ATAN(x)+0.5_r8*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-10.15_r8*ZoL)**r3
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN((1.0_r8+2.0_r8*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        bulk_psiu=(1.0_r8-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        bulk_psiu=-((1.0_r8+ZoL)+0.6667_r8*(ZoL-14.28_r8)/              &
     &            EXP(cff)+8.525_r8)
      END IF
      RETURN
      END FUNCTION bulk_psiu

      FUNCTION bulk_psit (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!                                                                      !
!
      USE mod_kinds
!
!  Function result
!
      real(r8) :: bulk_psit
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.5_r8
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-34.15_r8*ZoL)**r3
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN((1.0_r8+2.0_r8*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        bulk_psit=(1.0_r8-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        bulk_psit=-((1.0_r8+2.0_r8*ZoL)**1.5_r8+                        &
     &            0.6667_r8*(ZoL-14.28_r8)/EXP(cff)+8.525_r8)
      END IF

      RETURN
      END FUNCTION bulk_psit
#elif defined BULK_FLUXES2D && !defined CCSM_FLUXES2D
      implicit none
!
      PRIVATE
      PUBLIC  :: bulk_flux
!
      CONTAINS
!***********************************************************************
      SUBROUTINE bulk_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL bulk_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef ICE_MODEL
     &                     liold(ng),                                   &
     &                     linew(ng),                                   &
# endif
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % mask2,                            &
# endif
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif

      RETURN
      END SUBROUTINE bulk_flux
!
!***********************************************************************
      SUBROUTINE bulk_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liold, linew,                          &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
     &                           mask2,                                 &
# endif
     &                           Uwind, Vwind, sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
# ifdef ICE_MODEL
      integer, intent(in) :: liold
      integer, intent(in) :: linew
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: mask2(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: Uwind(LBi:,LBj:)
      real(r8), intent(out) :: Vwind(LBi:,LBj:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: mask2(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: Uwind(LBi:,LBj:)
      real(r8), intent(out) :: Vwind(LBi:,LBj:)
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, i, j, it, listp
      integer :: i0, j0
      integer, parameter :: IterMax = 3
      real(r8) :: cff
      real(r8), parameter :: tol = 0.001_r8
      real(r8), dimension(IminS:ImaxS) :: Wspeed
      real(r8), dimension(IminS:ImaxS) :: Cd
      real(r8), dimension(IminS:ImaxS) :: U10
      real(r8), dimension(IminS:ImaxS) :: U10o
      real(r8), dimension(IminS:ImaxS) :: diff
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy
      real(r8), parameter :: rho_air = 1.22         ! air density [kg/m^2]

# include "set_bounds.h"
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      cff = log(blk_ZW(ng)/10._r8)/vonKar
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          Wspeed(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          U10(i) = Wspeed(i)/(1 + cff*sqrt(1.15e-3_r8))
          U10o(i)=0.0_r8
          diff(i) = abs(U10(i)-U10o(i))
          Cd(i)=0.0_r8
        END DO
        DO while (maxval(diff(Istr-1:IendR)) > tol)
          DO i=Istr-1,IendR
            U10o(i) = u10(i)
            IF (U10o(i) .lt. 10.15) THEN
              Cd(i) = 1.15e-3_r8
            ELSE
              Cd(i) = 4.9e-4_r8 + 6.5e-5*U10o(i)
            END IF
            U10(i) = Wspeed(i)/(1+cff*sqrt(Cd(i)))
            diff(i) = abs(U10(i)-U10o(i))
          END DO
        END DO
        DO i=Istr-1,IendR
          Taux(i,j) = rho_air*Cd(i)*U10(i)*Uwind(i,j)
          Tauy(i,j) = rho_air*Cd(i)*U10(i)*Vwind(i,j)
        END DO
      END DO

      cff=0.5_r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
#undef MASK_HACK
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
#  ifdef MASK_HACK
! HACK for Bering strait test
          i0 = max(i-1,IstrR)
          sustr(i,j)= sustr(i,j) * max(mask2(i,j),mask2(i0,j))
#  endif
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
#  ifdef MASK_HACK
! HACK for Bering strait test
          j0 = max(j-1,JstrR)
          svstr(i,j)= svstr(i,j) * max(mask2(i,j),mask2(i,j0))
#  endif
# endif
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
      END SUBROUTINE bulk_flux_tile
#endif
      END MODULE bulk_flux_mod
