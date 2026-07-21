#include "cppdefs.h"
      MODULE bulk_flux_mod
#ifdef BULK_FLUXES
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
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
!      parameterization of the air-sea exchange of heat and            !
!      water vapor including the molecular constraints at              !
!      the interface, J. Atmos. Sci, 36, 1722-1735.                    !
!                                                                      !
!    Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea     !
!      surface roughness on the height and steepness of the waves.     !
!      J. Phys. Oceanogr., 31, 572-590.                                !
!                                                                      !
!    Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002: !
!      New evidence for a relation between wind stress and wave age    !
!      from measurements during ASGAMAGE. Bound.-Layer Meteor., 103,   !
!      409-438.                                                        !
!                                                                      !
!    Fairall, C.W., Bradley, E.F., Hare, J.E., Grachev, A.A. and       !
!      Edson, J.B., 2003. Bulk parameterization of air-sea fluxes:     !
!      Updates and verification for the COARE algorithm. Journal of    !
!      Climate, 16, 571-591.                                           !
!                                                                      !
!    Drennan, W. M., Graber, H. C., Hauser, D., & Quentin, C., 2003:   !
!      On the wave age dependence of wind stress over pure wind seas.  !
!      J. Geophys. Res: Oceans, 108(C3).                               !
!                                                                      !
!    Edson, J.B., Jampana, V., Weller, R.A., Bigorre, S.P.,            !
!      Plueddemann, A.J., Fairall, C.W., Miller, S.D., Mahrt, L.,      !
!      Vickers, D. and Hersbach, H., 2013. On the exchange of          !
!      momentum over the open ocean. Journal Physical Oceanography,    !
!      43, 1589-1610.                                                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
!  In this modified algorithm, the 'rain' variable is associated with  !
!  any phase of precipitation liquid or solid. For example, NARR data  !
!  includes rain, snow, freezing rain, etc.  So in most places where   !
!  'rain' was originally in this routine it has been replaced by       !
!  'ABS(rain)'. The negative sign/value is being used to indicate that !
!  the precipitation is of solid, rather than liquid form. The main    !
!  reason for distinguishing between liquid/solid precipitation is to  !
!  specify an additional latent heat flux to water when frozen water   !
!  falls upon it, or to ice when liquid water falls upon it. The       !
!  latent heat flux associated with the melting or freezing is         !
!  'hfus*rain'. It is only applied to the ocean surface when the       !
!  precipitation reaches the surface in solid form, it only applies to !
!  to the ice surface when the precipitation falls as a liquid (SMD).  !
!                                                                      !
# endif
!                                                                      !
!  Modified by Paul Goodman add EMINUSP equivalent salt flux (10/2004) !
!  Modified by Kate Hedstrom for COARE version 3.0 (03/2005).          !
!  Modified by Jim Edson to correct specific humidities.               !
!  Modified by John Wilkin to introduce WIND_MINUS_CURRENT (02/2018)   !
!  Modified by Fernando Pareja-Roman and John Wilkin to introduce      !
!      COARE 3.5 and correct the partition of vector stress (02/2024)  !
!      To recover older COAREv3.0 #define COARE_30                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      USE mod_ice
# endif
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
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
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17, __LINE__, MyFile)
# endif
      CALL bulk_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
     &                     liold(ng), linew(ng),                        &
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
# ifdef WIND_MINUS_CURRENT
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
# endif
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % Cd,                             &
# ifdef CLOUDS
     &                     FORCES(ng) % cloud,                          &
# endif
# if defined FRC_COUPLING && defined PBLH
     &                     FORCES(ng) % PBLh,                           &
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
# if defined ICE_MODEL  && defined ICE_BULK_FLUXES && \
     defined ICE_ALBEDO && defined SHORTWAVE
!!   &                     FORCES(ng) % albedo,                         &
     &                     FORCES(ng) % albedo_ice,                     &
# endif
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
     &                     ICE(ng) % Fi,                                &
     &                     ICE(ng) % Si,                                &
     &                     FORCES(ng) % sustr_ai,                       &
     &                     FORCES(ng) % svstr_ai,                       &
     &                     FORCES(ng) % sustr_ao,                       &
     &                     FORCES(ng) % svstr_ao,                       &
     &                     FORCES(ng) % Qnet_ao,                        &
#  ifdef ICE_THERMO
     &                     FORCES(ng) % Qnet_ai,                        &
     &                     FORCES(ng) % srflx_ice,                      &
     &                     FORCES(ng) % snow,                           &
#  endif
# endif
# ifdef EMINUSP
     &                     FORCES(ng) % evap,                           &
# endif
# ifdef WAV2OCN_FLUXES
     &                     FORCES(ng) % Tauocx,                         &
     &                     FORCES(ng) % Tauocy,                         &
# endif
     &                     FORCES(ng) % stflux,                         &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE bulk_flux
!
!***********************************************************************
      SUBROUTINE bulk_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
     &                           liold, linew,                          &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef WET_DRY
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
     &                           alpha, beta, rho, t,                   &
# ifdef WIND_MINUS_CURRENT
     &                           u, v,                                  &
# endif
     &                           Hair, Pair, Tair, Uwind, Vwind, Cd,    &
# ifdef CLOUDS
     &                           cloud,                                 &
# endif
# if defined FRC_COUPLING && defined PBLH
     &                           PBLh,                                  &
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
# if defined ICE_MODEL  && defined ICE_BULK_FLUXES && \
     defined ICE_ALBEDO && defined SHORTWAVE
     &                           albedo_ice,                            &
# endif
     &                           rain, lhflx, lrflx, shflx, srflx,      &
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
     &                           Fi, Si,                                &
     &                           sustr_ai, svstr_ai,                    &
     &                           sustr_ao, svstr_ao,                    &
     &                           Qnet_ao,                               &
#  ifdef ICE_THERMO
     &                           Qnet_ai,                               &
     &                           srflx_ice, snow,                       &
#  endif
# endif
# ifdef EMINUSP
     &                           evap,                                  &
# endif
# ifdef WAV2OCN_FLUXES
     &                           Tauocx, Tauocy,                        &
# endif
     &                           stflux, sustr, svstr)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      integer, intent(in) :: liold, linew
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in   ) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in   ) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in   ) :: alpha(LBi:,LBj:)
      real(r8), intent(in   ) :: beta(LBi:,LBj:)
      real(r8), intent(in   ) :: rho(LBi:,LBj:,:)
      real(r8), intent(in   ) :: t(LBi:,LBj:,:,:,:)
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in   ) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in   ) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in   ) :: Hair(LBi:,LBj:)
      real(r8), intent(in   ) :: Pair(LBi:,LBj:)
      real(r8), intent(in   ) :: Tair(LBi:,LBj:)
      real(r8), intent(in   ) :: Uwind(LBi:,LBj:)
      real(r8), intent(in   ) :: Vwind(LBi:,LBj:)
      real(r8), intent(inout) :: Cd(LBi:,LBj:)
#  ifdef CLOUDS
      real(r8), intent(in   ) :: cloud(LBi:,LBj:)
#  endif
#  if defined FRC_COUPLING && defined PBLH
      real(r8), intent(in   ) :: PBLh(LBi:,LBj:)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined DRENNAN
      real(r8), intent(in   ) :: Hwave(LBi:,LBj:)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN
      real(r8), intent(in   ) :: Pwave_top(LBi:,LBj:)
#  endif
#  if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
      real(r8), intent(in   ) :: Lwavep(LBi:,LBj:)
#  endif
#  if defined ICE_MODEL  && defined ICE_BULK_FLUXES && \
      defined ICE_ALBEDO && defined SHORTWAVE
      real(r8), intent(in   ) :: albedo_ice(LBi:,LBj:)
#  endif
      real(r8), intent(in   ) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflux(LBi:,LBj:,:)
#  if defined ICE_MODEL && defined ICE_BULK_FLUXES
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(out  ) :: sustr_ai(LBi:,LBj:)
      real(r8), intent(out  ) :: svstr_ai(LBi:,LBj:)
      real(r8), intent(out  ) :: sustr_ao(LBi:,LBj:)
      real(r8), intent(out  ) :: svstr_ao(LBi:,LBj:)
      real(r8), intent(out  ) :: Qnet_ao(LBi:,LBj:)
#   ifdef ICE_THERMO
      real(r8), intent(out  ) :: Qnet_ai(LBi:,LBj:)
      real(r8), intent(out  ) :: srflx_ice(LBi:,LBj:)
      real(r8), intent(out  ) :: snow(LBi:,LBj:)
#   endif
#  endif
#  ifdef EMINUSP
      real(r8), intent(out  ) :: evap(LBi:,LBj:)
#  endif
#  ifdef WAV2OCN_FLUXES
      real(r8), intent(in   ) :: Tauocx(LBi:,LBj:)
      real(r8), intent(in   ) :: Tauocy(LBi:,LBj:)
#  endif
      real(r8), intent(out  ) :: sustr(LBi:,LBj:)
      real(r8), intent(out  ) :: svstr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in   ) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in   ) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in   ) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in   ) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in   ) :: v(LBi:UBi,LBj:UBj,N(ng),3)
#  endif
      real(r8), intent(in   ) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Cd(LBi:UBi,LBj:UBj)
#  ifdef CLOUDS
      real(r8), intent(in   ) :: cloud(LBi:UBi,LBj:UBj)
#  endif
#  if defined FRC_COUPLING && defined PBLH
      real(r8), intent(in   ) :: PBLh(LBi:UBi,LBj:UBj)
#  endif
#  if defined COARE_TAYLOR_YELLAND || \
      defined DRENNAN
      real(r8), intent(in   ) :: Hwave(LBi:UBi,LBj:UBj)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN
      real(r8), intent(in   ) :: Pwave_top(LBi:UBi,LBj:UBj)
#  endif
#  if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST || \
      defined DRENNAN)
      real(r8), intent(in   ) :: Lwavep(LBi:UBi,LBj:UBj)
#  endif
#  if defined ICE_MODEL  && defined ICE_BULK_FLUXES && \
      defined ICE_ALBEDO && defined SHORTWAVE
      real(r8), intent(in   ) :: albedo_ice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in   ) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflux(LBi:UBi,LBj:UBj,NT(ng))
#  if defined ICE_MODEL && defined ICE_BULK_FLUXES
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(out  ) :: sustr_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: svstr_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: sustr_ao(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: svstr_ao(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: Qnet_ao(LBi:UBi,LBj:UBj)
#   ifdef ICE_THERMO
      real(r8), intent(out  ) :: Qnet_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: srflx_ice(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: snow(LBi:UBi,LBj:UBj)
#   endif
#  endif
#  ifdef EMINUSP
      real(r8), intent(out  ) :: evap(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WAV2OCN_FLUXES
      real(r8), intent(in   ) :: Tauocx(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: Tauocy(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out  ) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      integer :: li_stp
# endif
      integer, parameter :: IterMax = 3
!
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Hl, Hlw, Hscale, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      real(r8) :: Qsw_i, Qlw_i, Qlh_i, Qsh_i
      real(r8) :: Qsati, TiceK
      real(r8) :: Ce, Ch, dq_i, le_i, slp, vap_p_i
# endif
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: diffh, diffw, oL, upvel
      real(r8) :: twopi_inv, wet_bulb
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
# endif

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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uair
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vair
# ifdef ICE_THERMO
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Coef_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: RHS_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qai
# endif
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
# endif

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
# ifdef WIND_MINUS_CURRENT
!
!  Modify near-surface (2m or 10m) effective winds by subtracting the
!  ocean surface current (J Wilkin). See:
!
!  Bye, J.A.T. and J.-O. Wolff, 1999: Atmosphere-ocean momentum exchange
!     in general circulation models. J. Phys. Oceanogr. 29, 671-691.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,MIN(Iend+1,Lm(ng))
          Uair(i,j)=Uwind(i,j)-                                         &
     &              0.5_r8*(u(i,j,N(ng),nrhs)+u(i+1,j,N(ng),nrhs))
        END DO
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          Uair(Iend+1,j)=Uwind(Iend+1,j)-u(Iend,j,N(ng),nrhs)
        END IF
      END DO
      DO i=Istr-1,Iend+1
        DO j=Jstr-1,MIN(Jend+1,Mm(ng))
          Vair(i,j)=Vwind(i,j)-                                         &
     &              0.5_r8*(v(i,j,N(ng),nrhs)+v(i,j+1,N(ng),nrhs))
        END DO
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          Vair(i,Jend+1)=Vwind(i,Jend+1)-v(i,Jend,N(ng),nrhs)
        END IF
      END DO
# else
!
!  Load wind components to local arrays.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Uair(i,j)=Uwind(i,j)
          Vair(i,j)=Vwind(i,j)
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute Atmosphere-ocean fluxes using a bulk flux parameterization.
!-----------------------------------------------------------------------
!
      Hscale=rho0*Cp                              ! Celsius m/s to W/m2
      twopi_inv=0.5_r8/pi
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        li_stp=liold
      ELSE
        li_stp=linew
      END IF
# endif
!
      J_LOOP : DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
!  (At input the flux data is converted to Celsius m/s. Here, we need to
!   convert back to W/m2 using 'Hscale').
!
          Wmag(i)=SQRT(Uair(i,j)*Uair(i,j)+Vair(i,j)*Vair(i,j))+eps
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)                                 ! Celsius
          TairK(i)=TairC(i)+273.16_r8                        ! Kelvin
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)                   ! Celsius
          TseaK(i)=TseaC(i)+273.16_r8                        ! Kelvin
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
# if defined ICE_MODEL && defined ICE_BULK_FLUXES && defined ICE_THERMO
          srflx_ice(i,j)=SRad(i,j)
# endif
          Tcff(i)=alpha(i,j)
          Scff(i)=beta(i,j)
!
!  Initialize.
!
          delTc(i)=0.0_r8
          delQc(i)=0.0_r8
          LHeat(i,j)=lhflx(i,j)*Hscale     ! latent heat data, if any
          SHeat(i,j)=shflx(i,j)*Hscale     ! sensible heat data, if any
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
          cff=(0.7859_r8+0.03477_r8*TairC(i))/                          &
     &        (1.0_r8+0.00412_r8*TairC(i))
          e_sat=10.0_r8**cff   ! saturation vapor pressure (hPa or mbar)
          vap_p=e_sat*RH       ! water vapor pressure (hPa or mbar)
          cff2=TairK(i)*TairK(i)*TairK(i)
          cff1=cff2*TairK(i)
          LRad(i,j)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_r8-0.05_r8*SQRT(vap_p))*                &
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
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*                                            &
     &           (cff/(PairM-0.378_r8*cff+eps))        !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
!
!  Vapor Pressure reduced for salinity (Kraus and Businger, 1994, pp42).
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
# ifdef REGCM_COUPLING
          Hlv(i,j)=2.5104E+6_r8           ! to be consistent with RegCM
# else
          Hlv(i,j)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
# endif
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
          Cd(i,j)=(vonKar/LOG(blk_ZW(ng)/Zo10(i)))**2
!
!  Compute Richardson number.
# if defined FRC_COUPLING && defined PBLH
!  Use atmosphere boundary layer thickness from coupled model.
# endif
!
          Ct(i)=vonKar/LOG(blk_ZT(ng)/ZoT10(i))  ! T transfer coefficient
          CC(i)=vonKar*Ct(i)/Cd(i,j)
          delTc(i)=0.0_r8
# ifdef COOL_SKIN
          delTc(i)=blk_dter
# endif
# if defined FRC_COUPLING && defined PBLH
          Ribcu(i)=-blk_ZW(ng)/(PBLh(i,j)*0.004_r8*blk_beta**3)
# else
          Ribcu(i)=-blk_ZW(ng)/(blk_Zabl*0.004_r8*blk_beta**3)
# endif
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

# if defined COARE_30
!
!  Momentum roughness length based on COARE 3.0 (Fairall et al 2003, JPO)
!
          IF (delW(i).gt.18.0_r8) THEN
            charn(i)=0.018_r8
          ELSE IF ((10.0_r8.lt.delW(i)).and.(delW(i).le.18.0_r8)) THEN
            charn(i)=0.011_r8+                                          &
     &               0.125_r8*(0.018_r8-0.011_r8)*(delW(i)-10.0_r8)
          ELSE
            charn(i)=0.011_r8
          END IF

# else
!
! Momentum roughness length based on COARE 3.5 (Edson et al 2013, JPO)
!
          charn(i)=MIN(0.028_r8,-0.005_r8+0.0017_r8*delW(i))

# endif

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
# if defined DRAGLIM_DAVIS
            ZoW(i)=MAX(ZoW(i),1.27E-7_r8) ! Davis et al. (2008)
            ZoW(i)=MIN(ZoW(i),2.85E-3_r8)
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
#if defined COARE_30
!
! Moisture and thermal roughness lengths, COARE 3.0
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6_r8)
#else
!
! In COARE 3.5, moisture and thermal roughness lengths are adjusted
! to reflect COARE 3.0 fluxes (J. Edson, pers. comm, 11/2023)
!
            ZoQ(i)=MIN(1.6e-4_r8,5.8e-5_r8/(Rr(i)**0.72_r8))
# endif
            ZoT(i)=ZoQ(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoL(i)=vonKar*g*blk_ZW(ng)*                                 &
     &             (Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                     &
     &                        0.61_r8*TairK(i)*Qstar(i))/               &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &              (1.0_r8+0.61_r8*Q(i))+eps)
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
# if defined FRC_COUPLING && defined PBLH
!  Use atmosphere boundary layer thickness from coupled model.
# endif
!
            Bf=-g/(TairK(i)+0.61_r8*Qstar(i))*                          &
     &         Wstar(i)*(Tstar(i)+0.61_r8*TairK(i)*Qstar(i))
            IF (Bf.gt.0.0_r8) THEN
# if defined FRC_COUPLING && defined PBLH
              Wgus(i)=blk_beta*(Bf*PBLh(i,j))**r3
# else
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
# endif
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
            Clam=16.0_r8*g*blk_Cpw*(rhoSea(i)*blk_visw)**3.0_r8/        &
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
     &               (Clam*Qbouy/(Wstar(i)+eps)**4.0_r8)**0.75_r8)**r3
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
!  Compute nondimensional transfer coefficients for stress (Cd),
!  sensible heat (Ch), and latent heat (Ce).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd(i,j)=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
!         Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps) ! used prior to 02/2024
# if defined ICE_MODEL && defined ICE_BULK_FLUXES && defined ICE_THERMO
!         Ch=Wstar(i)*Tstar(i)/(-Wspeed*delT(i)+0.0098_r8*blk_ZT(ng))
!         Ce=Wstar(i)*Qstar(i)/(-Wspeed*delQ(i)+eps)
          Ch = 1.75E-3_r8                     ! Tale Bakken Ulfsby, 2025
          Ce = 1.75E-3_r8
# endif
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
     &          (rhoAir(i)*blk_Cpa+eps)
          cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
          wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/        &
     &                                     (blk_Cpa*diffh))
          Hsr=ABS(rain(i,j))*wet_bulb*blk_Cpw*                          &
     &        ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
# ifndef REGCM_COUPLING
          SHeat(i,j)=(Hs+Hsr)
# endif
# ifdef MASKING
          SHeat(i,j)=SHeat(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          SHeat(i,j)=SHeat(i,j)*rmask_wet(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
          Hl=-Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)

# if defined ICE_MODEL && defined ICE_BULK_FLUXES && defined ICE_THERMO
!
!  Subtract loss of heat associated with latent heat flux from ocean to
!  melting snow (SMD).
!
          IF (rain(i,j).lt.0.0_r8) THEN
            Hl=Hl+3.347E+5_r8*rain(i,j)
          END IF
# endif
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61_r8*Wstar(i)*Qstar(i)-                             &
     &          (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
          Hlw=rhoAir(i)*Hlv(i,j)*upvel*Q(i)
# ifndef REGCM_COUPLING
          LHeat(i,j)=(Hl+Hlw)
# endif
# ifdef MASKING
          LHeat(i,j)=LHeat(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          LHeat(i,j)=LHeat(i,j)*rmask_wet(i,j)
# endif
!
!  Horizontal momentum flux (N/m2) due to rain (kg/m2/s) impact
!  traveling at 85% of air velocity
!
          Taur=0.85_r8*ABS(rain(i,j))*Wmag(i)
!
!  Sum of stresses (N/m2) expressed as friction velocities.
!  Normalization by Wmag (m/s) so that stress components (Taux,
!  Tauy in N/m2) are cff times the respective wind components
!  (Uair,Vair in m/s).
!
          cff=rhoAir(i)*(Wstar(i)*Wstar(i)+Taur/rhoAir(i))/(Wmag(i)+eps)
!
!  In prior versions incorrect normalization included gustiness in
!  Wspeed causing the vector stress magnitude to be too be low by ~3%
!  (Corrected J. Wilkin and F. Parela-Roman 02/2024)
!
# ifdef WAV2OCN_FLUXES
          Taux(i,j)=Tauocx(i,j)
# else
          Taux(i,j)=cff*Uair(i,j)
# endif
# ifdef MASKING
          Taux(i,j)=Taux(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          Taux(i,j)=Taux(i,j)*rmask_wet(i,j)
# endif
# ifdef WAV2OCN_FLUXES
          Tauy(i,j)=Tauocy(i,j)
# else
          Tauy(i,j)=cff*Vair(i,j)
# endif
# ifdef MASKING
          Tauy(i,j)=Tauy(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          Tauy(i,j)=Tauy(i,j)*rmask_wet(i,j)
# endif
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
          Taux_Ice(i,j)=rhoAir(i)*Cd_ai(ng)*Uwind(i,j)*Wspeed
          Tauy_Ice(i,j)=rhoAir(i)*Cd_ai(ng)*Vwind(i,j)*Wspeed

#  ifdef ICE_THERMO
!
!  Correct ocean net short-wave radiation for ice cover and convert to
!  kinematic units.
!
          Coef_Ice(i,j)=0.0_r8
          RHS_Ice(i,j)=0.0_r8
!
!  Sensible heat flux over ice.
!
          IF (Si(i,j,li_stp,isAice).gt.min_ai(ng)) THEN
            TiceK=Fi(i,j,icIsst)+273.16_r8
          ELSE
            TiceK=t(i,j,N(ng),nrhs,itemp)+273.16_r8
          ENDIF
          cff=rhoAir(i)*blk_Cpa*Ch*delW(i)                   ! (W/m2)/K
          Qsh_i=cff*(TairK(i)+0.0098_r8*blk_ZT(ng)-TiceK)    ! W/m2
          Coef_Ice(i,j)=Coef_Ice(i,j)+cff                    ! (W/m2)/K
          RHS_Ice(i,j)=RHS_Ice(i,j)+                                    &
     &                 cff*(TairK(i)+0.0098_r8*blk_ZT(ng))   ! W/m2
!
!  Add in the sensible heat transfer associated with warm rain or cold
!  snow contacting the ice usinf the wet bulb factor in the sensible
!  heat flux for water calculation above. This, of course, is a bit off
!  if the actual surface is snow or meltwater covered. But in those
!  cases there are inaccuracies in other terms as well (SMD).
!
          cff=wet_bulb*2093.0_r8
          Qsh_i=Qsh_i+                                                  &
     &          ABS(rain(i,j))*cff*                                     &
     &          (TiceK-TairC(i)+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
          Coef_Ice(i,j)=Coef_Ice(i,j) -                                 &
     &                  ABS(rain(i,j))*cff
!
!  Use TairK here, unlike in above calculation of Hsr, since
!  Fi(:,:,icIsst) will be calculated in Kelvin (SMD, 10.20.20).
!
          RHS_Ice(i,j)=RHS_Ice(i,j)+                                    &
     &                 ABS(rain(i,j))*cff*                              &
     &                 (-TairK(i)+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
!
!  Latent heat flux over ice.
!
          le_i=2.834E+6_r8
!
!  Compute saturation specific humidity over the ice (Qsati) from
!  Parkinson and Washington (1979).
!
          slp=Pair(i,j)*100.0_r8           ! Convert back to Pa from mb
          cff=611.0_r8*                                                 &
     &        10.0_r8**(9.5_r8*(TiceK-273.16_r8)/(TiceK-7.66_r8))
!!        Qsati=eps*cff/(slp-(1.0_r8-.62197_r8)*cff)
          Qsati=0.62197_r8*cff/(slp-(1.0_r8-0.62197_r8)*cff)
          dq_i=Q(i)-Qsati
          Qlh_i=rhoAir(i)*Ce*delW(i)*                                   &
     &          le_i*dq_i
!
!  Add in the Latent heat gain from rain freezing on ice/snow surface.
!  If meltwater is already present then assume rain just pools,
!  otherwise rain may freeze or ice may melt (SMD).
!
          IF ((rain(i,j).gt.0.0_r8).and.                                &
     &        (Si(i,j,li_stp,isHmel).eq.0.0_r8)) THEN
            Qlh_i=Qlh_i+3.347E+5_r8*rain(i,j)
          END IF
          RHS_Ice(i,j)=RHS_Ice(i,j)+Qlh_i
!
!  Short-wave radiation (W/m**2) to ice.
!
          Qsw_i=(1.0_r8-albedo_ice(i,j))*srflx_ice(i,j)
!
!  A question remains as to what fraction of the solar radiation should
!  be considered to contribute to the surface temperature of the
!  ice/snow/meltwater?  Leaving unchanged for now, assuming entire
!  'Qsw_i' contributes to ice surface temperature although it should
!  probably be something less (SMD).
!
          RHS_Ice(i,j)=RHS_Ice(i,j)+Qsw_i
!
!  Net longwave radiation from ice.
!
          cff=Q(i)/(1.0_r8+Q(i))
          vap_p_i=slp*cff/(0.62197_r8+cff)

#   ifdef LONGWAVE
          cff1=0.39_r8-0.005_r8*SQRT(vap_p_i)
          cff2=1.0_r8-0.65_r8*cloud(i,j)
          cff3=StefBo*emmiss*(TairK(i)*TairK(i)*TairK(i))
          cff4=cff3*TairK(i)
          Qlw_i=cff4*cff1*cff2+                                         &
     &          4.0_r8*cff3*(TiceK-TairK(i))
          Coef_Ice(i,j)=Coef_Ice(i,j)+                                  &
     &                  4.0_r8*cff3
          RHS_Ice(i,j)=RHS_Ice(i,j)-                                    &
     &                 cff4*(cff1*cff2-4.0_r8)

#   elif defined LONGWAVE_OUT
!
!  Option where the incoming longwave radiation (downward) is given as
!  input (lrflx). The following should be consistent, with the implicit
!  (pseudo) approach to solving for 'Fi(:,:,icIsst). Here, 'Qlw_i' is
!  the flux from ice to atmosphere, thus the negative sign (SMD).
!
          cff3=StefBo*emmiss*(TiceK**3)
          cff4=cff3*TiceK
          Qlw_i=-(lrflx(i,j)*Hscale-cff4)
          Coef_Ice(i,j)=Coef_Ice(i,j)+                                  &
     &                  4.0_r8*cff3
          RHS_Ice(i,j)=RHS_Ice(i,j)+                                    &
     &                 lrflx(i,j)*Hscale+3.0_r8*cff4

#   else
          cff3=StefBo*emmiss*(TairK(i)**3)
          Qlw_i=lrflx(i,j)*Hscale
          Coef_Ice(i,j)=Coef_Ice(i,j)+                                  &
     &                  4.0_r8*cff3
          RHS_Ice(i,j)=RHS_Ice(i,j)-                                    &
     &                 lrflx(i,j)*Hscale-                               &
     &                 4.0_r8*cff3*(TiceK-2.0_r8*TairK(i))
#   endif
!
!  Net heat flux from ice to atmosphere. The heat equation terms are
!  computed with quantities in Kelvin to facilitate the computation of
!  seaice surface temperature, Fi(:,:,icIsst), in Kelivin in module
!  'ice_mk.h' and then converted to Celsius afterwards (SMD).
!
           Qai(i,j)=-Qsh_i-Qlh_i-Qsw_i+Qlw_i
#  endif
# endif
        END DO
      END DO J_LOOP
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
!  Note that stflux(:,:,isalt) is the E-P flux. The fresh water flux
!  is positive out of the ocean and the salt flux is positive into the
!  ocean. It is  multiplied by surface salinity when computing state
!  variable stflx(:,:,isalt) in "set_vbc.F".
!
      Hscale=1.0_r8/(rho0*Cp)                     ! W/2 to Celsius m/s
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
# if defined ALBEDO && defined SHORTWAVE
!!        srflx(i,j) = srflx(i,j)*(1.0_r8-albedo(i,j))
# endif
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
!
!  Limit shortwave radiation to be used in computing penetrative
!  radiation by presence of ice alone.
!
          cff1=1.0_r8/(Si(i,j,li_stp,isAice)+eps)
          ice_thick(i,j) =Si(i,j,li_stp,isHice)*cff1
          snow_thick(i,j)=Si(i,j,li_stp,isHsno)*cff1

          IF (ice_thick(i,j).le.0.1_r8) THEN
            cff2=ice_thick(i,j)*5.0_r8
          ELSE
            cff2=0.1_r8*(5.0_r8+(ice_thick(i,j)-0.1_r8))
          ENDIF

          cff2=cff2+snow_thick(i,j)*20.0_r8
          srflx(i,j)=(1.0_r8-Si(i,j,li_stp,isAice))*srflx(i,j)+         &
     &               Si(i,j,li_stp,isAice)*srflx(i,j)*EXP(-cff2)
!
!  Net precipitation (m/s). Negative values of 'rain' indicate frozen
!  precipitation (SMD).
!
          IF (rain(i,j).ge.0.0_r8) THEN
            snow(i,j)=0.0_r8
          ELSE
            snow(i,j)=ABS(rain(i,j))/SnowDryRho(ng)
          END IF
# endif
!
!  Load surface fluxes to state arrays.
!
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
          Fi(i,j,icQcon)=Coef_Ice(i,j)                       !(W/m2)/K
          Fi(i,j,icQrhs)=RHS_Ice(i,j)                        ! W/m2
# endif
          lrflx(i,j)=LRad(i,j)*Hscale                        ! C m/s
          lhflx(i,j)=-LHeat(i,j)*Hscale                      ! C m/s
          shflx(i,j)=-SHeat(i,j)*Hscale                      ! C m/s
          stflux(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                     &
     &                       lhflx(i,j)+shflx(i,j))          ! C m/s

# if defined ICE_MODEL && defined ICE_BULK_FLUXES
          Qnet_ao(i,j)=-stflux(i,j,itemp)/Hscale             ! W/m2
          Qnet_ai(i,j)=Qai(i,j)                              ! W/m2
# endif
# ifdef MASKING
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask_wet(i,j)
# endif
# ifdef EMINUSP
          evap(i,j)=LHeat(i,j)/(Hlv(i,j)+eps)                ! kg/m2/s
#  if defined ICE_MODEL && defined ICE_BULK_FLUXES
          evap(i,j)=(1.0_r8-Si(i,j,li_stp,isAice))*evap(i,j)
#  endif
#  ifdef MASKING
          evap(i,j)=evap(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          evap(i,j)=evap(i,j)*rmask_wet(i,j)
#  endif
#  ifdef SALINITY
          stflux(i,j,isalt)=cff*(evap(i,j)-rain(i,j))        ! m/s
#   ifdef MASKING
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask(i,j)
#   endif
#   ifdef WET_DRY
!#    if !defined EMINUSP_SHIMA
!!    Evap is masked on the dry cells above.
!!    We do not want rain to be masked on dry cells.
!          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask_wet(i,j)
!#    endif
!
#    if defined EMINUSP_SHIMA
!    try new way
!    if cell is dry, then eminusp can only be neg
          IF (rmask_wet(i,j).lt.0.5_r8) THEN
            IF (stflux(i,j,isalt).gt.0.0_r8) THEN
              stflux(i,j,isalt)=stflux(i,j,isalt)*rmask_wet(i,j)
            END IF
          END IF
#    else
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask_wet(i,j)
#    endif

#   endif
#  endif
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
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# endif
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
          sustr_ao(i,j)=sustr(i,j)
          sustr_ai(i,j)=Taux_Ice(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
          svstr_ao(i,j)=svstr(i,j)
          svstr_ai(i,j)=Tauy_Ice(i,j)
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
     &                          stflux(:,:,itemp))

# if defined ICE_MODEL && defined ICE_BULK_FLUXES
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Qnet_ao)

#  ifdef ICE_THERMO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Qnet_ai)

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icQcon))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icQrhs))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          snow)
#  endif

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr_ai)

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr_ai)
# endif

# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
#  ifdef SALINITY
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflux(:,:,isalt))
#  endif
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)

        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
# if defined ICE_MODEL && defined ICE_BULK_FLUXES
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                       &
     &                          sustr_ao)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                       &
     &                          svstr_ao)
# endif
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflux(:,:,itemp))

#  if defined ICE_MODEL && defined ICE_BULK_FLUXES
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, Qnet_ao, sustr_ai, svstr_ai)

#   ifdef ICE_THERMO
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
                          Fi(:,:,icQcon),                               &
     &                    Fi(:,:,icQrhs),                               &
     &                    Qnet_ai, snow)
#   endif
#  endif

#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap)
#   ifdef SALINITY
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflux(:,:,isalt))
#   endif
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)

# if defined ICE_MODEL && defined ICE_BULK_FLUXES
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr_ao, svstr_ao)
#  endif
# endif
!
      RETURN
      END SUBROUTINE bulk_flux_tile
!
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
      real(dp), intent(in) :: pi
      real(r8), intent(in) :: ZoL
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
!
      RETURN
      END FUNCTION bulk_psiu
!
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
      real(dp), intent(in) :: pi
      real(r8), intent(in) :: ZoL
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
!
      RETURN
      END FUNCTION bulk_psit
!
#endif
      END MODULE bulk_flux_mod
