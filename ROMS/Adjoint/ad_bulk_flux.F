#include "cppdefs.h"
      MODULE ad_bulk_flux_mod
#if defined ADJOINT && defined BULK_FLUXES_NOT_YET
!
!svn $Id: ad_bulk_flux.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
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
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_bulk_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bulk_flux (ng, tile)
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
      CALL wclock_on (ng, iADM, 17)
# endif
      CALL ad_bulk_flux_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS  JminS, JmaxS,               &
     &                        nrhs(ng),                                 &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
     &                        MIXING(ng) % alpha,                       &
     &                        MIXING(ng) % ad_alpha,                    &
     &                        MIXING(ng) % beta,                        &
     &                        MIXING(ng) % ad_beta,                     &
     &                        OCEAN(ng) % rho,                          &
     &                        OCEAN(ng) % ad_rho,                       &
     &                        OCEAN(ng) % t,                            &
     &                        OCEAN(ng) % ad_t,                         &
     &                        FORCES(ng) % Hair,                        &
     &                        FORCES(ng) % Pair,                        &
     &                        FORCES(ng) % Tair,                        &
     &                        FORCES(ng) % Uwind,                       &
     &                        FORCES(ng) % Vwind,                       &
# ifdef CLOUDS
     &                        FORCES(ng) % cloud,                       &
# endif
# ifdef BBL_MODEL_NOT_YET
     &                        FORCES(ng) % Awave,                       &
     &                        FORCES(ng) % Pwave,                       &
# endif
     &                        FORCES(ng) % rain,                        &
     &                        FORCES(ng) % lhflx,                       &
     &                        FORCES(ng) % ad_lhflx,                    &
     &                        FORCES(ng) % lrflx,                       &
     &                        FORCES(ng) % ad_lrflx,                    &
     &                        FORCES(ng) % shflx,                       &
     &                        FORCES(ng) % ad_shflx,                    &
     &                        FORCES(ng) % srflx,                       &
     &                        FORCES(ng) % stflx,                       &
     &                        FORCES(ng) % ad_stflx,                    &
# ifdef EMINUSP
     &                        FORCES(ng) % evap,                        &
     &                        FORCES(ng) % ad_evap,                     &
# endif
     &                        FORCES(ng) % ad_sustr,                    &
     &                        FORCES(ng) % ad_svstr)

# ifdef PROFILE
      CALL wclock_off (ng, iADM, 17)
# endif
      RETURN
      END SUBROUTINE ad_bulk_flux
!
!***********************************************************************
      SUBROUTINE ad_bulk_flux_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS          &
     &                              nrhs,                               &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
     &                              alpha, ad_alpha,                    &
     &                              beta, ad_beta,                      &
     &                              rho, ad_rho, t, ad_t,               &
     &                              Hair, Pair, Tair,                   &
     &                              Uwind, Vwind,                       &
# ifdef CLOUDS
     &                              cloud,                              &
# endif
# ifdef BBL_MODEL_NOT_YET
     &                              Awave, Pwave,                       &
# endif
     &                              rain,                               &
     &                              lhflx, ad_lhflx,                    &
     &                              lrflx, ad_lrflx,                    &
     &                              shflx, ad_shflx,                    &
     &                              srflx,                              &
     &                              stflx, ad_stflx,                    &
# ifdef EMINUSP
     &                              evap, ad_evap,                      &
# endif
     &                              ad_sustr, ad_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bulk_flux_mod, ONLY : bulk_psiu, bulk_psit
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
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
      real(r8), intent(inout) :: ad_alpha(LBi:,LBj:)
      real(r8), intent(inout) :: ad_beta(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
#  ifdef BBL_MODEL_NOT_YET
      real(r8), intent(in) :: Awave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: rain(LBi:,LBj:)

      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_shflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)

#  ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:,LBj:)
      real(r8), intent(inout) :: ad_evap(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: ad_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_beta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
#  ifdef BBL_MODEL_NOT_YET
      real(r8), intent(in) :: Awave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: ad_lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))

#  ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, IterA, i, j, k
      integer, parameter :: IterMax = 3

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Cd, Hl, Hlw, Hscale, Hs, Hsr, IER
      real(r8) :: ad_Bf, ad_Cd, ad_Hl, ad_Hlw, ad_Hsr, ad_Hs
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL

      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(r8) :: twopi_inv, wet_bulb
      real(r8) :: ad_wet_bulb, ad_Wspeed
      real(r8) :: fac, ad_fac, fac1, ad_fac1, fac2, ad_fac2
      real(r8) :: fac3, ad_fac3, ad_cff, ad_upvel
      real(r8) :: adfac, adfac1, adfac2
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
      real(r8) :: ad_Clam,  ad_Hcool, ad_Hsb, ad_Hlb
      real(r8) :: ad_Qbouy, ad_Qcool, ad_lambd
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
      real(r8), dimension(IminS:ImaxS) :: delW1
      real(r8), dimension(IminS:ImaxS) :: L
      real(r8), dimension(IminS:ImaxS) :: L10
      real(r8), dimension(IminS:ImaxS) :: Lwave
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

      real(r8), dimension(IminS:ImaxS) :: Wstar1
      real(r8), dimension(IminS:ImaxS) :: Tstar1
      real(r8), dimension(IminS:ImaxS) :: Qstar1
      real(r8), dimension(IminS:ImaxS) :: delTc1

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy

      real(r8), dimension(IminS:ImaxS) :: ad_CC
      real(r8), dimension(IminS:ImaxS) :: ad_charn
      real(r8), dimension(IminS:ImaxS) :: ad_Ct
      real(r8), dimension(IminS:ImaxS) :: ad_Cd10
      real(r8), dimension(IminS:ImaxS) :: ad_Ct10
      real(r8), dimension(IminS:ImaxS) :: ad_Cwet
      real(r8), dimension(IminS:ImaxS) :: ad_delQ
      real(r8), dimension(IminS:ImaxS) :: ad_delQc
      real(r8), dimension(IminS:ImaxS) :: ad_delT
      real(r8), dimension(IminS:ImaxS) :: ad_delTc
      real(r8), dimension(IminS:ImaxS) :: ad_delW
      real(r8), dimension(IminS:ImaxS) :: ad_L
      real(r8), dimension(IminS:ImaxS) :: ad_L10
      real(r8), dimension(IminS:ImaxS) :: ad_Q
      real(r8), dimension(IminS:ImaxS) :: ad_Qpsi
      real(r8), dimension(IminS:ImaxS) :: ad_Qsea
      real(r8), dimension(IminS:ImaxS) :: ad_Qstar
      real(r8), dimension(IminS:ImaxS) :: ad_rhoSea
      real(r8), dimension(IminS:ImaxS) :: ad_Ri
      real(r8), dimension(IminS:ImaxS) :: ad_Rr
      real(r8), dimension(IminS:ImaxS) :: ad_Scff
      real(r8), dimension(IminS:ImaxS) :: ad_Tcff
      real(r8), dimension(IminS:ImaxS) :: ad_Tpsi
      real(r8), dimension(IminS:ImaxS) :: ad_TseaC
      real(r8), dimension(IminS:ImaxS) :: ad_TseaK
      real(r8), dimension(IminS:ImaxS) :: ad_Tstar
      real(r8), dimension(IminS:ImaxS) :: ad_u10
      real(r8), dimension(IminS:ImaxS) :: ad_Wgus
      real(r8), dimension(IminS:ImaxS) :: ad_Wpsi
      real(r8), dimension(IminS:ImaxS) :: ad_Wstar
      real(r8), dimension(IminS:ImaxS) :: ad_Zetu
      real(r8), dimension(IminS:ImaxS) :: ad_Zo10
      real(r8), dimension(IminS:ImaxS) :: ad_ZoT10
      real(r8), dimension(IminS:ImaxS) :: ad_ZoL
      real(r8), dimension(IminS:ImaxS) :: ad_ZoQ
      real(r8), dimension(IminS:ImaxS) :: ad_ZoT
      real(r8), dimension(IminS:ImaxS) :: ad_ZoW
      real(r8), dimension(IminS:ImaxS) :: ad_ZWoL

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Tauy

# include "set_bounds.h"
!
!-------------------------------------------------------------------------
!  Initialize adjoint private variables.
!-------------------------------------------------------------------------
!
      ad_wet_bulb=0.0_r8
      ad_Wspeed=0.0_r8
      ad_Bf=0.0_r8
      ad_Cd=0.0_r8
      ad_Hl=0.0_r8
      ad_Hlw=0.0_r8
      ad_Hsr=0.0_r8
      ad_Hs=0.0_r8
      ad_fac=0.0_r8
      ad_fac1=0.0_r8
      ad_fac2=0.0_r8
      ad_fac3=0.0_r8
      ad_cff=0.0_r8
      ad_upvel=0.0_r8
# ifdef COOL_SKIN
      ad_Clam=0.0_r8
      ad_Hcool=0.0_r8
      ad_Hsb=0.0_r8
      ad_Hlb=0.0_r8
      ad_Qbouy=0.0_r8
      ad_Qcool=0.0_r8
      ad_lambd=0.0_r8
# endif
      DO i=IminS,ImaxS
        ad_CC(i)=0.0_r8
        ad_charn(i)=0.0_r8
        ad_Ct(i)=0.0_r8
        ad_Cd10(i)=0.0_r8
        ad_Ct10(i)=0.0_r8
        ad_Cwet(i)=0.0_r8
        ad_delQ(i)=0.0_r8
        ad_delQc(i)=0.0_r8
        ad_delT(i)=0.0_r8
        ad_delTc(i)=0.0_r8
        ad_delW(i)=0.0_r8
        ad_L(i)=0.0_r8
        ad_L10(i)=0.0_r8
        ad_Q(i)=0.0_r8
        ad_Qpsi(i)=0.0_r8
        ad_Qsea(i)=0.0_r8
        ad_Qstar(i)=0.0_r8
        ad_rhoSea(i)=0.0_r8
        ad_Ri(i)=0.0_r8
        ad_Rr(i)=0.0_r8
        ad_Scff(i)=0.0_r8
        ad_Tcff(i)=0.0_r8
        ad_Tpsi(i)=0.0_r8
        ad_TseaC(i)=0.0_r8
        ad_TseaK(i)=0.0_r8
        ad_Tstar(i)=0.0_r8
        ad_u10(i)=0.0_r8
        ad_Wgus(i)=0.0_r8
        ad_Wpsi(i)=0.0_r8
        ad_Wstar(i)=0.0_r8
        ad_Zetu(i)=0.0_r8
        ad_Zo10(i)=0.0_r8
        ad_ZoT10(i)=0.0_r8
        ad_ZoL(i)=0.0_r8
        ad_ZoQ(i)=0.0_r8
        ad_ZoT(i)=0.0_r8
        ad_ZoW(i)=0.0_r8
        ad_ZWoL(i)=0.0_r8
      END DO
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_Hlv(i,j)=0.0_r8
          ad_LHeat(i,j)=0.0_r8
          ad_LRad(i,j)=0.0_r8
          ad_SHeat(i,j)=0.0_r8
          ad_Taux(i,j)=0.0_r8
          ad_Tauy(i,j)=0.0_r8
        END DO
      END DO
!
      twopi_inv=0.5_r8/pi
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_sustr, tl_svstr)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_sustr, ad_svstr)
#  ifdef EMINUSP
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_evap,                                      &
!>   &                    tl_stflx(:,:,isalt))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_evap,                                   &
     &                       ad_stflx(:,:,isalt))
#  endif
!>    CALL mp_exchange2d (ng, tile, iTLM, 4,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_lrflx, tl_lhflx, tl_shflx,                 &
!>   &                    tl_stflx(:,:,itemp))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 4,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_lrflx, ad_lhflx, ad_shflx,              &
     &                       ad_stflx(:,:,itemp))
# endif
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_svstr)
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_svstr)
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_sustr)
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_sustr)
# ifdef EMINUSP
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_stflx(:,:,isalt))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_stflx(:,:,isalt))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_evap)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_evap)
# endif
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_stflx(:,:,itemp))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_stflx(:,:,itemp))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_shflx)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_shflx)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_lhflx)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_lhflx)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_lrflx)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_lrflx)
      END IF
!
!=======================================================================
!  Compute adjoint surface net heat flux and surface wind stress.
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
!  to (psu m/s) for stflx(isalt) in "set_vbc.F".
!


!
!  Compute adjoint of kinematic, surface wind stress (m2/s2).
!
      cff=0.5_r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
# ifdef MASKING
!>        tl_sustr(i,j)=tl_sustr(i,j)*umask(i,j)
          ad_sustr(i,j)=ad_sustr(i,j)*umask(i,j)
# endif
!>        tl_sustr(i,j)=cff*(tl_Taux(i-1,j)+tl_Taux(i,j))
          adfac=cff*ad_sustr(i,j)
          ad_Taux(i-1,j)=ad_Taux(i-1,j)+adfac
          ad_Taux(i  ,j)=ad_Taux(i  ,j)+adfac
# if !(defined AD_SENSITIVITY      || defined ADJUST_WSTRESS   || \
       defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
       defined SO_SEMI)
          ad_sustr(i,j)=0.0_r8
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
# ifdef MASKING
!>        tl_svstr(i,j)=tl_svstr(i,j)*vmask(i,j)
          ad_svstr(i,j)=ad_svstr(i,j)*vmask(i,j)
# endif
!>        tl_svstr(i,j)=cff*(tl_Tauy(i,j-1)+tl_Tauy(i,j))
          adfac=cff*ad_svstr(i,j)
          ad_Tauy(i,j-1)=ad_Tauy(i,j-1)+adfac
          ad_Tauy(i,j  )=ad_Tauy(i,j  )+adfac
# if !(defined AD_SENSITIVITY      || defined ADJUST_WSTRESS   || \
       defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
       defined SO_SEMI)
          ad_svstr(i,j)=0.0_r8
# endif
        END DO
      END DO

!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          Hlv(i,j)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
        END DO
      END DO

      Hscale=1.0_r8/(rho0*Cp)
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
# ifdef MASKING
#   ifdef EMINUSP
!>        tl_evap(i,j)=tl_evap(i,j)*rmask(i,j)
          ad_evap(i,j)=ad_evap(i,j)*rmask(i,j)
!>        tl_stflx(i,j,isalt)=tl_stflx(i,j,isalt)*rmask(i,j)
          ad_stflx(i,j,isalt)=ad_stflx(i,j,isalt)*rmask(i,j)
#   endif
!>        tl_stflx(i,j,itemp)=tl_stflx(i,j,itemp)*rmask(i,j)
          ad_stflx(i,j,itemp)=ad_stflx(i,j,itemp)*rmask(i,j)
# endif
# ifdef EMINUSP
!>        tl_stflx(i,j,isalt)=cff*tl_evap(i,j)
          ad_evap(i,j)=ad_evap(i,j)+cff*ad_stflx(i,j,isalt)
#  if !(defined AD_SENSITIVITY      || defined ADJUST_STFLUX    || \
        defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
        defined SO_SEMI)
          ad_stflx(i,j,isalt)=0.0_r8
#  endif
!>        tl_evap(i,j)=tl_LHeat(i,j)/Hlv(i,j)-                          &
!>   &                               tl_Hlv(i,j)*evap(i,j)/Hlv(i,j)
          adfac=ad_evap(i,j)/Hlv(i,j)
          ad_LHeat(i,j)=ad_LHeat(i,j)+adfac
          ad_Hlv(i,j)=ad_Hlv(i,j)-adfac*evap(i,j)
#  if !(defined AD_SENSITIVITY   || defind IS4DVAR_SENSITIVITY || \
        defined OPT_OBSERVATIONS || defined SO_SEMI)
          ad_evap(i,j)=0.0_r8
#  endif
# endif
!>        tl_stflx(i,j,itemp)=(tl_lrflx(i,j)+                           &
!>   &                      tl_lhflx(i,j)+tl_shflx(i,j))
          ad_lrflx(i,j)=ad_lrflx(i,j)+ad_stflx(i,j,itemp)
          ad_lhflx(i,j)=ad_lhflx(i,j)+ad_stflx(i,j,itemp)
          ad_shflx(i,j)=ad_shflx(i,j)+ad_stflx(i,j,itemp)
# if !(defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
       defined OPT_OBSERVATIONS || defined SO_SEMI)
          ad_stflx(i,j,itemp)=0.0_r8
# endif
!>        tl_shflx(i,j)=-tl_SHeat(i,j)*Hscale
          ad_SHeat(i,j)=ad_SHeat(i,j)-ad_shflx(i,j)*Hscale
# if !(defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
       defined OPT_OBSERVATIONS || defined SO_SEMI)
          ad_shflx(i,j)=0.0_r8
# endif
!>        tl_lhflx(i,j)=-tl_LHeat(i,j)*Hscale
          ad_LHeat(i,j)=ad_LHeat(i,j)-ad_lhflx(i,j)*Hscale
# if !(defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
       defined OPT_OBSERVATIONS || defined SO_SEMI)
          ad_lhflx(i,j)=0.0_r8
# endif
!>        tl_lrflx(i,j)=tl_LRad(i,j)*Hscale
          ad_LRad(i,j)=ad_LRad(i,j)+ad_lrflx(i,j)*Hscale
          ad_lrflx(i,j)=0.0_r8

        END DO
      END DO
!
!=======================================================================
!  Adjoint Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      DO j=Jstr-1,JendR
!
!  Compute Basic State quantities. required
!
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
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
!  1.0 at poles to 0.5 at the equator).
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
          Qair(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!  Compute specific humidity, Q (kg/kg).
!
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &            EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
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
     &          (TairK(i)*delW(i)*delW(i))
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
# ifdef BBL_MODEL_NOT_YET
          Cwave(i)=g*Pwave(i,j)*twopi_inv
          Lwave(i)=Cwave(i)*Pwave(i,j)
# endif
        END DO
!
!  Iterate until convergence. It usually converges within 3 iterations.
!
        DO Iter=1,IterMax
          DO i=Istr-1,IendR
# ifdef BBL_MODEL_NOT_YET
!
!  Use wave info if we have it, two different options.
!
#  ifdef WIND_WAVES
            ZoW(i)=(25._r8/pi)*Lwave(i)*(Wstar(i)/Cwave(i))**4.5+       &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  else
            ZoW(i)=1200._r8*Awave(i,j)*(Awave(i,j)/Lwave(i))**4.5+      &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  endif
# else
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6)
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
     &              (1.0_r8+(Clam*Qbouy/(Wstar(i)+eps)**4)**0.75_r8)**r3
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

!
!-----------------------------------------------------------------------
!  Compute Adjoint of Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
        DO i=Istr-1,IendR
!
!  Compute transfer coefficients for momentum (Cd).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
            Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
            diffw=2.11E-5_r8*(TairK(i)/273.16_r8)**1.94_r8
            diffh=0.02411_r8*(1.0_r8+TairC(i)*                          &
     &                        (3.309E-3_r8-1.44E-6_r8*TairC(i)))/       &
     &            (rhoAir(i)*blk_Cpa)
            cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
            wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/      &
     &                                       (blk_Cpa*diffh))
            Hsr=rain(i,j)*wet_bulb*blk_Cpw*                             &
     &          ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
            SHeat(i,j)=(Hs+Hsr)
# ifdef MASKING
            SHeat(i,j)=SHeat(i,j)*rmask(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
            Hl=-Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
            upvel=-1.61_r8*Wstar(i)*Qstar(i)-                           &
     &            (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
            Hlw=rhoAir(i)*Hlv(i,j)*upvel*Q(i)
            LHeat(i,j)=(Hl+Hlw)
# ifdef MASKING
            LHeat(i,j)=LHeat(i,j)*rmask(i,j)
# endif
!
!  Adjoint wind stress components (N/m2), Tau.
!

# ifdef MASKING
!>          tl_Tauy(i,j)=tl_Tauy(i,j)*rmask(i,j)
            ad_Tauy(i,j)=ad_Tauy(i,j)*rmask(i,j)
# endif
!>          tl_Tauy(i,j)=tl_cff*Vwind(i,j)
            ad_cff=ad_cff+ad_Tauy(i,j)*Vwind(i,j)
            ad_Tauy(i,j)=0.0_r8
# ifdef MASKING
!>          tl_Taux(i,j)=tl_Taux(i,j)*rmask(i,j)
            ad_Taux(i,j)=ad_Taux(i,j)*rmask(i,j)
# endif
!>          tl_Taux(i,j)=tl_cff*Uwind(i,j)
            ad_cff=ad_cff+ad_Taux(i,j)*Uwind(i,j)
            ad_Taux(i,j)=0.0_r8

!>          tl_cff=rhoAir(i)*(tl_Cd*Wspeed+Cd*tl_Wspeed)
            adfac=rhoAir(i)*ad_cff
            ad_Cd=ad_Cd+Wspeed*adfac
            ad_Wspeed=ad_Wspeed+Cd*adfac
            ad_cff=0.0_r8
!
!  Adjoint turbulent latent heat flux (W/m2), Hl.
!


# ifdef MASKING
!>          tl_LHeat(i,j)=tl_LHeat(i,j)*rmask(i,j)
            ad_LHeat(i,j)=ad_LHeat(i,j)*rmask(i,j)
# endif
!>          tl_LHeat(i,j)=(tl_Hl+tl_Hlw)
            ad_Hl=ad_Hl+ad_LHeat(i,j)
            ad_Hlw=ad_Hlw+ad_LHeat(i,j)
            ad_LHeat(i,j)=0.0_r8

!>          tl_Hlw=rhoAir(i)*Q(i)*(tl_Hlv(i,j)*upvel+                   &
!>   &                         Hlv(i,j)*tl_upvel)
            adfac=rhoAir(i)*Q(i)*ad_Hlw
            ad_Hlv(i,j)=ad_Hlv(i,j)+upvel*adfac
            ad_upvel=ad_upvel+Hlv(i,j)*adfac
            ad_Hlw=0.0_r8

!>          tl_upvel=-1.61_r8*                                          &
!>   &            (tl_Wstar(i)*Qstar(i)+Wstar(i)*tl_Qstar(i))-          &
!>   &            (1.0_r8+1.61_r8*Q(i))*(tl_Wstar(i)*Tstar(i)+          &
!>   &                                Wstar(i)*tl_Tstar(i))/TairK(i)
            adfac=-1.61_r8*ad_upvel
            adfac1=-(1.0_r8+1.61_r8*Q(i))*ad_upvel/TairK(i)
            ad_Wstar(i)=ad_Wstar(i)+Qstar(i)*adfac+Tstar(i)*adfac1
            ad_Qstar(i)=ad_Qstar(i)+Wstar(i)*adfac
            ad_Tstar(i)=ad_Tstar(i)+Wstar(i)*adfac1
            ad_upvel=0.0_r8

!>          tl_Hl=-tl_Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)-             &
!>   &             Hlv(i,j)*rhoAir(i)*(tl_Wstar(i)*Qstar(i)+            &
!>   &                                     Wstar(i)*tl_Qstar(i) )
            adfac=Hlv(i,j)*rhoAir(i)*ad_Hl
            ad_Hlv(i,j)=ad_Hlv(i,j)-rhoAir(i)*Wstar(i)*Qstar(i)*ad_Hl
            ad_Wstar(i)=ad_Wstar(i)-Qstar(i)*adfac
            ad_Qstar(i)=ad_Qstar(i)-Wstar(i)*adfac
            ad_Hl=0.0_r8
!
!  Adjoint sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!


# ifdef MASKING
!>          tl_SHeat(i,j)=tl_SHeat(i,j)*rmask(i,j)
            ad_SHeat(i,j)=ad_SHeat(i,j)*rmask(i,j)
# endif

!>          tl_SHeat(i,j)=(tl_Hs+tl_Hsr)
            ad_Hs=ad_Hs+ad_SHeat(i,j)
            ad_Hsr=ad_Hsr+ad_SHeat(i,j)
            ad_SHeat(i,j)=0.0_r8

!>          tl_Hsr=Hsr*tl_wet_bulb/wet_bulb+                            &
!>   &          rain(i,j)*wet_bulb*blk_Cpw*                             &
!>   &          (tl_TseaC(i)+tl_Qsea(i)*Hlv(i,j)/blk_Cpa+               &
!>   &          (Qsea(i)-Q(i))*tl_Hlv(i,j)/blk_Cpa)
            adfac=rain(i,j)*wet_bulb*blk_Cpw*ad_Hsr
            adfac1=adfac/blk_Cpa
            ad_wet_bulb=ad_wet_bulb+Hsr*ad_Hsr/wet_bulb
            ad_TseaC(i)=ad_TseaC(i)+adfac
            ad_Qsea(i)=ad_Qsea(i)+Hlv(i,j)*adfac1
            ad_Hlv(i,j)=ad_Hlv(i,j)+(Qsea(i)-Q(i))*adfac1
            ad_Hsr=0.0_r8

!>          tl_wet_bulb=-tl_fac*wet_bulb*wet_bulb
            ad_fac=-wet_bulb*wet_bulb*ad_wet_bulb
            ad_wet_bulb=0.0_r8

!>          tl_fac=0.622_r8*diffw*(tl_cff*Hlv(i,j)+cff*tl_Hlv(i,j))/    &
!>   &                                      (blk_Cpa*diffh)
            adfac=0.622_r8*diffw*ad_fac/(blk_Cpa*diffh)
            ad_cff=ad_cff+Hlv(i,j)*adfac
            ad_Hlv(i,j)=ad_Hlv(i,j)+cff*adfac
            ad_fac=0.0_r8

!>          tl_cff=Qair(i)*tl_Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
            ad_Hlv(i,j)=ad_Hlv(i,j)+Qair(i)*ad_cff/                     &
     &                                     (blk_Rgas*TairK(i)*TairK(i))
            ad_cff=0.0_r8
!
!  Adjoint turbulent sensible heat flux (W/m2), Hs.
!

!>          tl_Hs=-blk_Cpa*rhoAir(i)*(tl_Wstar(i)*Tstar(i)+             &
!>   &                                Wstar(i)*tl_Tstar(i))
            adfac=-blk_Cpa*rhoAir(i)*ad_Hs
            ad_Wstar(i)=ad_Wstar(i)+Tstar(i)*adfac
            ad_Tstar(i)=ad_Tstar(i)+Wstar(i)*adfac
            ad_Hs=0.0_r8

!
!  Adjoint transfer coefficients for momentum (Cd).
!
!>        tl_Cd=2.0_r8*(tl_Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)        &
!>   &          -tl_Wspeed*Wspeed*Cd/(Wspeed*Wspeed+eps))
          adfac=2.0_r8*ad_Cd/(Wspeed*Wspeed+eps)
          ad_Wstar(i)=ad_Wstar(i)+Wstar(i)*adfac
          ad_Wspeed=ad_Wspeed-Wspeed*Cd*adfac
          ad_Cd=0.0_r8

          IF (Wspeed.ne.0.0_r8) THEN
!>          tl_Wspeed=tl_Wgus(i)*Wgus(i)/Wspeed
            ad_Wgus(i)=ad_Wgus(i)+Wgus(i)*ad_Wspeed/Wspeed
            ad_Wspeed=0.0_r8
          ELSE
!>          tl_Wspeed=0.0_r8
            ad_Wspeed=0.0_r8
          END IF

        END DO
!
!  Adjoint Iteration until convergence. It usually converges within
!  3 iterations.
!
        DO Iter=IterMax,1,-1
!
!  Recompute appropriate basic state variables again prior to
!  each adjoint iteration.
!
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
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
!  1.0 at poles to 0.5 at the equator).
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
          Qair(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!  Compute specific humidity, Q (kg/kg).
!
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &            EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
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
     &          (TairK(i)*delW(i)*delW(i))
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
# ifdef BBL_MODEL_NOT_YET
          Cwave(i)=g*Pwave(i,j)*twopi_inv
          Lwave(i)=Cwave(i)*Pwave(i,j)
# endif
        END DO
!
!  Iterate here to obtain the appropriate basic state variables.
!
        DO IterA=1,Iter-1
          DO i=Istr-1,IendR
# ifdef BBL_MODEL_NOT_YET
!
!  Use wave info if we have it, two different options.
!
#  ifdef WIND_WAVES
            ZoW(i)=(25._r8/pi)*Lwave(i)*(Wstar(i)/Cwave(i))**4.5+       &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  else
            ZoW(i)=1200._r8*Awave(i,j)*(Awave(i,j)/Lwave(i))**4.5+      &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  endif
# else
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6)
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
     &              (1.0_r8+(Clam*Qbouy/(Wstar(i)+eps)**4)**0.75_r8)**r3
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

          DO i=Istr-1,IendR
!
!   Now complete the computation of the appropriate basic state quantities
!   for this iteration.
!
# ifdef BBL_MODEL_NOT_YET
!
!  Use wave info if we have it, two different options.
!
#  ifdef WIND_WAVES
            ZoW(i)=(25._r8/pi)*Lwave(i)*(Wstar(i)/Cwave(i))**4.5+       &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  else
            ZoW(i)=1200._r8*Awave(i,j)*(Awave(i,j)/Lwave(i))**4.5+      &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
#  endif
# else
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6)
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
            delQc(i)=Cwet(i)*delTc(i)
# endif
!
!  Compute wind scaling parameters, Wstar.
!
            Wstar1(i)=MAX(eps,delW(i)*vonKar/                           &
     &               (LOG(blk_ZW(ng)/ZoW(i))-Wpsi(i)))
            Tstar1(i)=-(delT(i)-delTc(i))*vonKar/                       &
     &               (LOG(blk_ZT(ng)/ZoT(i))-Tpsi(i))
            Qstar1(i)=-(delQ(i)-delQc(i))*vonKar/                       &
     &               (LOG(blk_ZQ(ng)/ZoQ(i))-Qpsi(i))
!
!  Compute gustiness in wind speed.
!
            Bf=-g/TairK(i)*                                             &
     &         Wstar1(i)*(Tstar1(i)+0.61_r8*TairK(i)*Qstar1(i))
            IF (Bf.gt.0.0_r8) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
            ELSE
              Wgus(i)=0.2_r8
            END IF
            delW1(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
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
            Hsb=-rhoAir(i)*blk_Cpa*Wstar1(i)*Tstar1(i)
            Hlb=-rhoAir(i)*Hlv(i,j)*Wstar1(i)*Qstar1(i)
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
     &             (1.0_r8+(Clam*Qbouy/(Wstar1(i)+eps)**4)**0.75_r8)**r3
              Hcool=lambd*blk_visw/(SQRT(rhoAir(i)/rhoSea(i))*          &
     &                              Wstar1(i)+eps)
              delTc1(i)=Qcool*Hcool/blk_tcw
            ELSE
              delTc1(i)=0.0_r8
            END IF
            delQc(i)=Cwet(i)*delTc1(i)
!
!-----------------------------------------------------------------------
!  Adjoint Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Adjoint temperature and moisture change.
!


!>          tl_delQc(i)=tl_Cwet(i)*delTc1(i)+Cwet(i)*tl_delTc(i)
            ad_Cwet(i)=ad_Cwet(i)+delTc1(i)*ad_delQc(i)
            ad_delTc(i)=ad_delTc(i)+Cwet(i)*ad_delQc(i)
            ad_delQc(i)=0.0_-r8

            IF ((Qcool.gt.0.0_r8).and.(Qbouy.gt.0.0_r8)) THEN
              fac1=SQRT(rhoAir(i)/rhoSea(i))
              fac2=fac1*Wstar1(i)+eps
              Hcool=lambd*blk_visw/fac2

!>            tl_delTc(i)=(tl_Qcool*Hcool+Qcool*tl_Hcool)/blk_tcw
              adfac=ad_delTc(i)/blk_tcw
              ad_Qcool=ad_Qcool+Hcool*adfac
              ad_Hcool=ad_Hcool+Qcool*adfac
              ad_delTc(i)=0.0_r8

!>            tl_Hcool=tl_lambd*blk_visw/fac2-tl_fac2*Hcool/fac2
              adfac=ad_Hcool/fac2
              ad_lambd=ad_lambd+blk_visw*adfac
              ad_fac2=-Hcool*adfac
              ad_Hcool=0.0_r8

              IF (fac1.ne.0.0_r8) THEN
!>                tl_fac1=-0.5_r8*tl_rhoSea(i)*fac1/rhoSea(i)
                  ad_rhoSea(i)=ad_rhoSea(i)-0.5*fac1*ad_fac1/rhoSea(i)
                  ad_fac1=0.0_r8
              ELSE
!>                tl_fac1=0.0_r8
                  ad_fac1=0.0_r8
              END IF

              fac1=(Wstar1(i)+eps)**4
              fac2=Clam*Qbouy
              fac3=(fac2/fac1)**0.75_r8
              lambd=6.0_r8/(1.0_r8+fac3)**r3

!>            tl_lambd=-r3*6.0_r8*tl_fac3/(1.0_r8+fac3)**(r3+1.0_r8)
              ad_fac3=-r3*6.0_r8*ad_lambd/(1.0_r8+fac3)**(r3+1.0_r8)
              ad_lambd=0.0_r8

!>            tl_fac3=0.75_r8*(fac2/fac1)**(-0.25_r8)*                  &
!>   &              (tl_fac2/fac1-tl_fac1*fac2/(fac1*fac1))
              adfac=0.75_r8*(fac2/fac1)**(-0.25_r8)*ad_fac3
              ad_fac2=adfac/fac1
              ad_fac1=-fac2*adfac/(fac1*fac1)
              ad_fac3=0.0_r8

!>            tl_fac2=tl_Clam*Qbouy+Clam*tl_Qbouy

              ad_Clam=ad_Clam+Qbouy*ad_fac2
              ad_Qbouy=ad_Qbouy+Clam*ad_fac2
              ad_fac2=0.0_r8

!>            tl_fac1=4.0_r8*tl_Wstar(i)*(Wstar1(i)+eps)**3

              ad_Wstar(i)=ad_Wstar(i)+4.0_r8*ad_fac1*(Wstar1(i)+eps)**3
              ad_fac1=0.0_r8

            ELSE
!>            tl_delTc(i)=0.0_r8
              ad_delTc(i)=0.0_r8
            END IF
!
!  Adjoint of Total cooling at the interface.
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
            Hsb=-rhoAir(i)*blk_Cpa*Wstar1(i)*Tstar1(i)
            Hlb=-rhoAir(i)*Hlv(i,j)*Wstar1(i)*Qstar1(i)
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

!>          tl_Qbouy=tl_Tcff(i)*Qcool+Tcff(i)*tl_Qcool+                 &
!>   &                  (tl_Scff(i)*Hlb*blk_Cpw+                        &
!>   &                  Scff(i)*tl_Hlb*blk_Cpw)/Hlv(i,j)-               &
!>   &               tl_Hlv(i,j)*Scff(i)*Hlb*blk_Cpw/(Hlv(i,j)*Hlv(i,j))

            adfac=ad_Qbouy*blk_Cpw/Hlv(i,j)
            ad_Tcff(i)=ad_Tcff(i)+Qcool*ad_Qbouy
            ad_Qcool=ad_Qcool+Tcff(i)*ad_Qbouy
            ad_Scff(i)=ad_Scff(i)+Hlb*adfac
            ad_Hlb=ad_Hlb+Scff(i)*adfac
            ad_Hlv(i,j)=ad_Hlv(i,j)-                                    &
     &                 Scff(i)*Hlb*adfac/Hlv(i,j)
            ad_Qbouy=0.0_r8


!>          tl_Qcool=tl_LRad(i,j)+tl_Hsb+tl_Hlb
            ad_LRad(i,j)=ad_LRad(i,j)+ad_Qcool
            ad_Hsb=ad_Hsb+ad_Qcool
            ad_Hlb=ad_Hlb+ad_Qcool
            ad_Qcool=0.0_r8
!
! Adjoint  Backgound sensible and latent heat.
!
!>          tl_Hlb=-rhoAir(i)*(tl_Hlv(i,j)*Wstar1(i)*Qstar1(i)+         &
!>   &                         Hlv(i,j)*(tl_Wstar(i)*Qstar1(i)+         &
!>   &                                 Wstar1(i)*tl_Qstar(i))

            adfac=-rhoAir(i)*ad_Hlb
            adfac1=Hlv(i,j)*adfac
            ad_Hlv(i,j)=ad_Hlv(i,j)+Wstar1(i)*Qstar1(i)*adfac
            ad_Wstar(i)=ad_Wstar(i)+Qstar1(i)*adfac1
            ad_Qstar(i)=ad_Qstar(i)+Wstar1(i)*adfac1
            ad_Hlb=0.0_r8

!>          tl_Hsb=-rhoAir(i)*blk_Cpa*(tl_Wstar(i)*Tstar1(i)+           &
!>   &                                 Wstar1(i)*tl_Tstar(i))
            adfac=-rhoAir(i)*blk_Cpa*ad_Hsb
            ad_Wstar(i)=ad_Wstar(i)+Tstar1(i)*adfac
            ad_Tstar(i)=ad_Tstar(i)+Wstar1(i)*adfac
            ad_Hsb=0.0_r8

!
! Adjoint  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
!>            tl_Clam=48.0_r8*g*blk_Cpw*tl_rhoSea(i)*blk_visw*          &
!>   &                (rhoSea(i)*blk_visw)**2/                          &
!>   &                (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))

              ad_rhoSea(i)=ad_rhoSea(i)+                                &
     &                     48.0_r8*g*blk_Cpw*ad_Clam*blk_visw*          &
     &                     (rhoSea(i)*blk_visw)**2/                     &
     &                     (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))
              ad_Clam=0.0_r8
# endif

!
! Adjoint  Compute gustiness in wind speed.
!
            IF (delW1(i).ne.0.0_r8)THEN
!>             tl_delW(i)=tl_Wgus(i)*Wgus(i)/delW1(i)
               ad_Wgus(i)=ad_Wgus(i)+Wgus(i)*ad_delW(i)/delW1(i)
               ad_delW(i)=0.0_r8
            ELSE
!>             tl_delW(i)=0.0_r8
               ad_delW(i)=0.0_r8
            END IF

            IF (Bf.gt.0.0_r8) THEN
!>            tl_Wgus(i)=r3*blk_beta*tl_Bf*blk_Zabl*                    &
!>   &                                      (Bf*blk_Zabl)**(r3-1.0_r8)
              ad_Bf=ad_Bf+r3*blk_beta*ad_Wgus(i)*blk_Zabl*              &
     &                                      (Bf*blk_Zabl)**(r3-1.0_r8)
              ad_Wgus(i)=0.0_r8
            ELSE
!>            tl_Wgus(i)=0.0_r8
              ad_Wgus(i)=0.0_r8
            END IF

!>          tl_Bf=-g/TairK(i)*(                                         &
!>   &         tl_Wstar(i)*(Tstar1(i)+0.61_r8*TairK(i)*Qstar1(i))+      &
!>   &         Wstar1(i)*(tl_Tstar(i)+0.61_r8*TairK(i)*tl_Qstar(i)) )
            adfac=-g/TairK(i)*ad_Bf
            adfac1=Wstar1(i)*adfac
            ad_Wstar(i)=ad_Wstar(i)+                                    &
     &                   (Tstar1(i)+0.61_r8*TairK(i)*Qstar1(i))*adfac
            ad_Tstar(i)=ad_Tstar(i)+adfac1
            ad_Qstar(i)=ad_Qstar(i)+0.61_r8*TairK(i)*adfac1
            ad_Bf=0.0_r8

!
! Adjoint of wind scaling parameters, Wstar.
!

            fac1=blk_ZQ(ng)/ZoQ(i)
            fac2=LOG(fac1)

!>          tl_Qstar(i)=-(tl_delQ(i)-tl_delQc(i))*vonKar/(fac2-Qpsi(i))-&
!>   &                    (tl_fac2-tl_Qpsi(i))*Qstar1(i)/(fac2-Qpsi(i))
            adfac=-ad_Qstar(i)*vonKar/(fac2-Qpsi(i))
            adfac1=-ad_Qstar(i)*Qstar1(i)/(fac2-Qpsi(i))
            ad_delQ(i)=ad_delQ(i)+adfac
            ad_delQc(i)=ad_delQc(i)-adfac
            ad_fac2=adfac1
            ad_Qpsi(i)=ad_Qpsi(i)-adfac1
            ad_Qstar(i)=0.0_r8

!>          tl_fac2=tl_fac1/fac1
            ad_fac1=ad_fac2/fac1
            ad_fac2=0.0_r8

!>          tl_fac1=-tl_ZoQ(i)*fac1/ZoQ(i)
            ad_ZoQ(i)=ad_ZoQ(i)-fac1*ad_fac1/ZoQ(i)
            ad_fac1=0.0_r8

            fac1=blk_ZT(ng)/ZoT(i)
            fac2=LOG(fac1)

!>          tl_Tstar(i)=-(tl_delT(i)-tl_delTc(i))*vonKar/(fac2-Tpsi(i))-&
!>   &           (tl_fac2-tl_Tpsi(i))*Tstar1(i)/(fac2-Tpsi(i))
            adfac=-ad_Tstar(i)*vonKar/(fac2-Tpsi(i))
            adfac1=-ad_Tstar(i)*Tstar1(i)/(fac2-Tpsi(i))
            ad_delT(i)=ad_delT(i)+adfac
            ad_delTc(i)=ad_delTc(i)-adfac
            ad_fac2=adfac1
            ad_Tpsi(i)=ad_Tpsi(i)-adfac1
            ad_Tstar(i)=0.0_r8

!>          tl_fac2=tl_fac1/fac1
            ad_fac1=ad_fac1+ad_fac2/fac1
            ad_fac2=0.0_r8

!>          tl_fac1=-tl_ZoT(i)*fac1/ZoT(i)
            ad_ZoT(i)=ad_ZoT(i)-fac1*ad_fac1/ZoT(i)
            ad_fac1=0.0_r8

            fac1=blk_ZW(ng)/ZoW(i)
            fac2=LOG(fac1)
            fac3=delW(i)*vonKar/(fac2-Wpsi(i))

!>          tl_Wstar(i)=(0.5_r8-SIGN(0.5_r8,eps-fac3))*tl_fac3
            ad_fac3=ad_Wstar(i)*(0.5_r8-SIGN(0.5_r8,eps-fac3))
            ad_Wstar(i)=0.0_r8

!>          tl_fac3=tl_delW(i)*vonKar/(fac2-Wpsi(i))-                   &
!>   &             (tl_fac2-tl_Wpsi(i))*fac3/(fac2-Wpsi(i))
            adfac=ad_fac3/(fac2-Wpsi(i))
            adfac1=fac3*adfac
            ad_delW(i)=ad_delW(i)+vonKar*adfac
            ad_fac2=-adfac1
            ad_Wpsi(i)=ad_Wpsi(i)+adfac1
            ad_fac3=0.0_r8

!>          tl_fac2=tl_fac1/fac1
            ad_fac1=ad_fac1+ad_fac2/fac1
            ad_fac2=0.0_r8

!>          tl_fac1=-tl_ZoW(i)*fac1/ZoW(i)
            ad_ZoW(i)=ad_ZoW(i)-fac1*ad_fac1/ZoW(i)
            ad_fac1=0.0_r8

# ifdef COOL_SKIN
!>          tl_delQc(i)=tl_Cwet(i)*delTc(i)+Cwet(i)*tl_delTc(i)
            ad_Cwet(i)=ad_Cwet(i)+delTc(i)*ad_delQc(i)
            ad_delTc(i)=ad_delTc(i)+Cwet(i)*ad_delQc(i)
            ad_delQc(i)=0.0_r8

!>          tl_Cwet(i)=0.622_r8*(tl_Hlv(i,j)*Qsea(i)+                   &
!>   &                                           Hlv(i,j)*tl_Qsea(i))/  &
!>   &              (blk_Rgas*TseaK(i)*TseaK(i))-                       &
!>   &               2.0_r8*blk_Rgas*tl_TseaK(i)*TseaK(i)*Cwet(i)/      &
!>   &              (blk_Rgas*TseaK(i)*TseaK(i))

            adfac=ad_Cwet(i)/(blk_Rgas*TseaK(i)*TseaK(i))
            adfac1=2.0_r8*blk_Rgas*TseaK(i)*Cwet(i)*adfac
            adfac2=0.622_r8*adfac
            ad_Hlv(i,j)=ad_Hlv(i,j)+Qsea(i)*adfac2
            ad_Qsea(i)=ad_Qsea(i)+Hlv(i,j)*adfac2
            ad_TseaK(i)=ad_TseaK(i)-adfac1
            ad_Cwet(i)=0.0_r8

# endif

!
!  Adjoint stability functions at Z/L.
!
            fac=blk_ZQ(ng)/L(i)
!>          tl_Qpsi(i)=tl_bulk_psit(tl_fac,fac,pi)
            ad_fac=ad_bulk_psit(ad_Qpsi(i),fac,pi)
            ad_Qpsi(i)=0.0_r8

!>          tl_fac=-tl_L(i)*fac/L(i)
            ad_L(i)=ad_L(i)-ad_fac*fac/L(i)
            ad_fac=0.0_r8

            fac=blk_ZT(ng)/L(i)
!>          tl_Tpsi(i)=tl_bulk_psit(tl_fac,fac,pi)
            ad_fac=ad_bulk_psit(ad_Tpsi(i),fac,pi)
            ad_Tpsi(i)=0.0_r8

!>          tl_fac=-tl_L(i)*fac/L(i)
            ad_L(i)=ad_L(i)-ad_fac*fac/L(i)
            ad_fac=0.0_r8

!>          tl_Wpsi(i)=tl_bulk_psiu(tl_ZoL(i),ZoL(i),pi)
            ad_ZoL(i)=ad_ZoL(i)+ad_bulk_psiu(ad_Wpsi(i),ZoL(i),pi)
            ad_Wpsi(i)=0.0_r8

!
!  Adjoint Monin-Obukhov stability parameter, Z/L.
!
!>          tl_L(i)=-tl_ZoL(i)*L(i)/(ZoL(i)+eps)
            ad_ZoL(i)=ad_ZoL(i)-L(i)*ad_L(i)/(ZoL(i)+eps)
            ad_L(i)=0.0_r8

!>          tl_ZoL(i)=vonKar*g*blk_ZW(ng)*                              &
!>   &             (tl_Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                  &
!>   &                        0.61_r8*TairK(i)*tl_Qstar(i))/            &
!>   &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
!>   &             (1.0_r8+0.61_r8*Q(i))+eps)-                          &
!>   &             2.0_r8*TairK(i)*tl_Wstar(i)*Wstar(i)*                &
!>   &             (1.0_r8+0.61_r8*Q(i))*ZoL(i)/                        &
!>   &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
!>   &             (1.0_r8+0.61_r8*Q(i))+eps) )
            adfac=vonKar*g*blk_ZW(ng)*ad_ZoL(i)/                        &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
            ad_Tstar(i)=ad_Tstar(i)+(1.0_r8+0.61_r8*Q(i))*adfac
            ad_Qstar(i)=ad_Qstar(i)+0.61_r8*TairK(i)*adfac
            ad_Wstar(i)=ad_Wstar(i)-2.0_r8*TairK(i)*Wstar(i)*           &
     &             (1.0_r8+0.61_r8*Q(i))*ZoL(i)*ad_ZoL(i)/              &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
            ad_ZoL(i)=0.0_r8

!>          tl_ZoT(i)=tl_ZoQ(i)
            ad_ZoQ(i)=ad_ZoQ(i)+ad_ZoT(i)
            ad_ZoT(i)=0.0_r8

!>          tl_ZoQ(i)=                                                  &
!>   &         -(0.5_r8-SIGN(0.5_r8,5.5e-5_r8/Rr(i)**0.6-1.15e-4_r8))   &
!>   &          *0.6_r8*5.5e-5_r8*tl_Rr(i)/Rr(i)**1.6
            ad_Rr(i)=ad_Rr(i)-                                          &
     &          (0.5_r8-SIGN(0.5_r8,5.5e-5_r8/Rr(i)**0.6-1.15e-4_r8))   &
     &          *0.6_r8*5.5e-5_r8*ad_ZoQ(i)/Rr(i)**1.6
            ad_ZoQ(i)=0.0_r8

!>          tl_Rr(i)=(tl_ZoW(i)*Wstar(i)+ZoW(i)*tl_Wstar(i))/VisAir(i)
            adfac=ad_Rr(i)/VisAir(i)
            ad_ZoW(i)=ad_ZoW(i)+Wstar(i)*adfac
            ad_Wstar(i)=ad_Wstar(i)+ZoW(i)*adfac
            ad_Rr(i)=0.0_r8

# ifdef BBL_MODEL_NOT_YET
#  ifdef WIND_WAVES
!>          tl_ZoW(i)=(25._r8/pi)*Lwave(i)*4.5_r8*tl_Wstar(i)*          &
!>   &                (Wstar(i)/Cwave(i))**3.5/Cwave(i)-                &
!>   &               tl_Wstar(i)*0.11_r8*VisAir(i)/                     &
!>   &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
            ad_Wstar(i)=ad_Wstar(i)+(25._r8/pi)*Lwave(i)*4.5_r8*        &
     &             (Wstar(i)/Cwave(i))**3.5*ad_ZoW(i)/Cwave(i)-         &
     &             0.11_r8*VisAir(i)*ad_ZoW(i)/                         &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
            ad_ZoW(i)0.0_r8
#  else
!>          tl_ZoW(i)=-tl_Wstar(i)*0.11_r8*VisAir(i)/                   &
!>   &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
            ad_Wstar(i)=ad_Wstar(i)-0.11_r8*VisAir(i)*ad_ZoW(i)/        &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
            ad_ZoW(i)0.0_r8
#  endif
# else
!>          tl_ZoW(i)=2.0_r8*charn(i)*tl_Wstar(i)*Wstar(i)/g-           &
!>   &             tl_Wstar(i)*0.11_r8*VisAir(i)/                       &
!>   &                           ((Wstar(i)+eps)*(Wstar(i)+eps))
            adfac=ad_ZoW(i)/g
            ad_Wstar(i)=ad_Wstar(i)+2.0_r8*charn(i)*Wstar(i)*adfac-     &
     &              0.11_r8*VisAir(i)*ad_ZoW(i)/                        &
     &                           ((Wstar(i)+eps)*(Wstar(i)+eps))
            ad_ZoW(i)=0.0_r8

# endif
          END DO
        END DO
!
!   End of Adjoint Iteration Loop
!

        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
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
!  1.0 at poles to 0.5 at the equator).
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
          Qair(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!  Compute specific humidity, Q (kg/kg).
!
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &            EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
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
     &          (TairK(i)*delW(i)*delW(i))
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


!
!  Adjoint of Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!

          IF (delW(i).gt.18.0_r8) THEN
!>          tl_charn(i)=0.0_r8
            ad_charn(i)=0.0_r8
          ELSE IF ((10.0_r8.lt.delW(i)).and.(delW(i).le.18.0_r8)) THEN
!>          tl_charn(i)=0.0_r8
            ad_charn(i)=0.0_r8
          ELSE
!>          tl_charn(i)=0.0_r8
            ad_charn(i)=0.0_r8
          END IF

!
!  Adjoint First guesses for Monon-Obukhov similarity scales.
!

          fac1=blk_ZQ(ng)/L10(i)
          fac2=LOG(blk_ZQ(ng)/ZoT10(i))

!>        tl_Qstar(i)=-(tl_delQ(i)-tl_delQc(i))*vonKar/                 &
!>   &             (fac2-bulk_psit(fac1,pi))-                           &
!>   &             (tl_fac2-tl_bulk_psit(tl_fac1,fac1,pi))*Qstar(i)/    &
!>   &             (fac2-bulk_psit(fac1,pi))

          adfac=ad_Qstar(i)*vonKar/                                     &
     &             (fac2-bulk_psit(fac1,pi))
          cff=Qstar(i)/(fac2-bulk_psit(fac1,pi))
          ad_delQ(i)=ad_delQ(i)-adfac
          ad_delQc(i)=ad_delQc(i)+adfac
          ad_fac2=-cff*ad_Qstar(i)
!>          ad_fac1=ad_bulk_psit(ad_Qstar(i),fac1,pi)*cff
          ad_fac1=ad_bulk_psit(cff*ad_Qstar(i),fac1,pi)
          ad_Qstar(i)=0.0_r8

!>        tl_fac2=-tl_ZoT10(i)/ZoT10(i)
          ad_ZoT10(i)=ad_ZoT10(i)-ad_fac2/ZoT10(i)
          ad_fac2=0.0_r8

!>        tl_fac1=-tl_L10(i)*fac1/L10(i)
          ad_L10(i)=ad_L10(i)-ad_fac1*fac1/L10(i)
          ad_fac1=0.0_r8

          fac1=blk_ZT(ng)/L10(i)
          fac2=LOG(blk_ZT(ng)/ZoT10(i))

!>        tl_Tstar(i)=-(tl_delT(i)-tl_delTc(i))*vonKar/                 &
!>   &             (fac2-bulk_psit(fac1,pi))-                           &
!>   &             (tl_fac2-tl_bulk_psit(tl_fac1,fac1,pi))*Tstar(i)/    &
!>   &             (fac2-bulk_psit(fac1,pi))
          adfac=ad_Tstar(i)*vonKar/                                     &
     &             (fac2-bulk_psit(fac1,pi))
          cff=Tstar(i)/(fac2-bulk_psit(fac1,pi))
          ad_delT(i)=ad_delT(i)-adfac
          ad_delTc(i)=ad_delTc(i)+adfac
          ad_fac2=-cff*ad_Tstar(i)
!>          ad_fac1=ad_bulk_psit(ad_Tstar(i),fac1,pi)*cff
          ad_fac1=ad_bulk_psit(cff*ad_Tstar(i),fac1,pi)
          ad_Tstar(i)=0.0_r8

!>        tl_fac2=-tl_ZoT10(i)/ZoT10(i)
          ad_ZoT10(i)=ad_ZoT10(i)-ad_fac2/ZoT10(i)
          ad_fac2=0.0_r8

!>        tl_fac1=-tl_L10(i)*fac1/L10(i)
          ad_L10(i)=ad_L10(i)-fac1*ad_fac1/L10(i)
          ad_fac1=0.0_r8

          fac1=blk_ZW(ng)/L10(i)
          fac2=LOG(blk_ZW(ng)/Zo10(i))

!>        tl_Wstar(i)=-(tl_fac2-tl_bulk_psiu(tl_fac1,fac1,pi))*Wstar(i) &
!>   &                                 /(fac2-bulk_psiu(fac1,pi))
          cff=Wstar(i)/(fac2-bulk_psiu(fac1,pi))
          ad_fac2=-cff*ad_Wstar(i)
!>          ad_fac1=cff*ad_bulk_psiu(ad_Wstar(i),fac1,pi)
          ad_fac1=ad_bulk_psiu(cff*ad_Wstar(i),fac1,pi)
          ad_Wstar(i)=0.0_r8

!>        tl_fac2=-tl_Zo10(i)/Zo10(i)
          ad_Zo10(i)=ad_Zo10(i)-ad_fac2/Zo10(i)
          ad_fac2=0.0_r8

!>        tl_fac1=-tl_L10(i)*fac1/L10(i)
          ad_L10(i)=ad_L10(i)-ad_fac1*fac1/L10(i)
          ad_fac1=0.0_r8

!
!  Adjoint Richardson number.
!

!>        tl_L10(i)=-L10(i)*L10(i)*tl_Zetu(i)/blk_ZW(ng)
          ad_Zetu(i)=ad_Zetu(i)-L10(i)*L10(i)*ad_L10(i)/blk_ZW(ng)
          ad_L10(i)=0.0_r8

          IF (Ri(i).lt.0.0_r8) THEN
!>          tl_Zetu(i)=(tl_CC(i)*Ri(i)+CC(i)*tl_Ri(i))/                 &
!>   &                                 (1.0_r8+Ri(i)/Ribcu(i))-         &
!>   &              (tl_Ri(i)/Ribcu(i))*Zetu(i)/(1.0_r8+Ri(i)/Ribcu(i))
            adfac=ad_Zetu(i)/(1.0_r8+Ri(i)/Ribcu(i))
            ad_CC(i)=ad_CC(i)+Ri(i)*adfac
            ad_Ri(i)=ad_Ri(i)+CC(i)*adfac-Zetu(i)*adfac/Ribcu(i)
            ad_Zetu(i)=0.0_r8
          ELSE
            fac=3.0_r8*Ri(i)/CC(i)
!>          tl_Zetu(i)=(tl_CC(i)*Ri(i)+CC(i)*tl_Ri(i))/(1.0_r8+fac)     &
!>   &                -tl_fac*Zetu(i)/(1.0_r8+fac)
            adfac=ad_Zetu(i)/(1.0_r8+fac)
            ad_CC(i)=ad_CC(i)+Ri(i)*adfac
            ad_Ri(i)=ad_Ri(i)+CC(i)*adfac
            ad_fac=-Zetu(i)*adfac
            ad_Zetu(i)=0.0_r8

!>          tl_fac=3.0_r8*tl_Ri(i)/CC(i)-tl_CC(i)*fac/CC(i)
            ad_Ri(i)=ad_Ri(i)+3.0_r8*ad_fac/CC(i)
            ad_CC(i)=ad_CC(i)-ad_fac*fac/CC(i)
            ad_fac=0.0_r8
          END IF

          fac=1.0/(TairK(i)*delW(i)*delW(i))

!>        tl_Ri(i)=-fac*g*blk_ZW(ng)*((tl_delT(i)-tl_delTc(i))+         &
!>   &                                 0.61_r8*TairK(i)*tl_delQ(i))
!>
          adfac=-fac*g*blk_ZW(ng)*ad_Ri(i)
          ad_delT(i)=ad_delT(i)+adfac
          ad_delTc(i)=ad_delTc(i)-adfac
          ad_delQ(i)=ad_delQ(i)+0.61_r8*TairK(i)*adfac
          ad_Ri(i)=0.0_r8

# ifdef COOL_SKIN
!>        tl_delTc(i)=0.0_r8
          ad_delTc(i)=0.0_r8
# endif

!>        tl_delTc(i)=0.0_r8
          ad_delTc(i)=0.0_r8

!>        tl_CC(i)=vonKar*tl_Ct(i)/Cd-tl_Cd*CC(i)/Cd
          adfac=ad_CC(i)/Cd
          ad_Ct(i)=ad_Ct(i)+vonKar*adfac
          ad_Cd=ad_Cd-CC(i)*adfac
          ad_CC(i)=0.0_r8

          fac=LOG(blk_ZT(ng)/ZoT10(i))

!>        tl_Ct(i)=-tl_fac*Ct(i)/fac
          ad_fac=-ad_Ct(i)*Ct(i)/fac

!>        tl_fac=-tl_ZoT10(i)/ZoT10(i)
          ad_ZoT10(i)=ad_ZoT10(i)-ad_fac/ZoT10(i)
          ad_fac=0.0_r8

!
!  Adjoint Neutral coefficients.
!
          fac=LOG(blk_ZW(ng)/Zo10(i))

!>        tl_Cd=-2.0_r8*tl_fac*Cd/fac
          ad_fac=-2.0_r8*ad_Cd/fac
          ad_Cd=0.0_r8

!>        tl_fac=-tl_Zo10(i)/Zo10(i)
          ad_Zo10(i)=ad_Zo10(i)-ad_fac/Zo10(i)
          ad_fac=0.0_r8

          fac=vonKar/Ct10(i)

!>        tl_ZoT10(i)=-tl_fac*ZoT10(i)
          ad_fac=-ad_ZoT10(i)*ZoT10(i)
          ad_fac=0.0_r8

!>        tl_fac=-tl_Ct10(i)*fac/Ct10(i)
          ad_Ct10(i)=ad_Ct10(i)-fac*ad_fac/Ct10(i)
          ad_fac=0.0_r8

!>        tl_Ct10(i)=-0.5_r8*tl_Cd10(i)*Ct10(i)/Cd10(i)
          ad_Cd10(i)=ad_Cd10(i)-0.5_r8*Ct10(i)*ad_Ct10(i)/Cd10(i)
          ad_Ct10(i)=0.0_r8

          fac=LOG(10.0_r8/Zo10(i))

!>        tl_Cd10(i)=-2.0_r8*tl_fac*Cd10(i)/fac
          ad_fac=-2.0_r8*Cd10(i)*ad_Cd10(i)/fac
          ad_Cd10(i)=0.0_r8

!>        tl_fac=-tl_Zo10(i)/Zo10(i)
          ad_Zo10(i)=ad_Zo10(i)-ad_fac/Zo10(i)
          ad_fac=0.0_r8

!>        tl_Zo10(i)=0.022_r8*tl_Wstar(i)*Wstar(i)/g-                   &
!>   &            tl_Wstar(i)*0.11_r8*VisAir(i)/(Wstar(i)*Wstar(i))

          ad_Wstar(i)=ad_Wstar(i)+0.022_r8*Wstar(i)*ad_Zo10(i)/g-       &
     &            0.11_r8*VisAir(i)*ad_Zo10(i)/(Wstar(i)*Wstar(i))
          ad_Zo10(i)=0.0_r8

!>        tl_Wstar(i)=0.035_r8*tl_u10(i)
          ad_u10(i)=ad_u10(i)+0.035_r8*ad_Wstar(i)
          ad_Wstar(i)=0.0_r8

!>        tl_u10(i)=0.0_r8

          ad_u10(i)=0.0_r8

!>        tl_ZoW(i)=0.0_r8
          ad_ZoW(i)=0.0_r8


!
!  Adjoint wind is measured relative to sea surface and include
!  gustiness.
!

!>        tl_delT(i)=tl_TseaC(i)
          ad_TseaC(i)=ad_TseaC(i)+ad_delT(i)
          ad_delT(i)=0.0_r8

!>        tl_delQ(i)=tl_Qsea(i)
          ad_Qsea(i)=ad_Qsea(i)+ad_delQ(i)
          ad_delQ(i)=0.0_r8

!>        tl_delW(i)=0.0_r8
          ad_delW(i)=0.0_r8

!>        tl_Wgus(i)=0.0_r8
          ad_Wgus(i)=0.0_r8

!
!  Adjoint latent heat of vaporization (J/kg) at sea surface, Hlv.
!
!>        tl_Hlv(i,j)=-0.00237_r8*tl_TseaC(i)*1.0E+6_r8
!>
          ad_TseaC(i)=ad_TseaC(i)-0.00237_r8*1.0E+6_r8*ad_Hlv(i,j)
          ad_Hlv(i,j)=0.0_r8
!
!  Adjoint Qsea (kg/kg) from vapor pressure.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
          cff=cff*0.98_r8
!>        tl_Qsea(i)=tl_fac*(1.0_r8+0.378_r8*cff/((PairM-0.378_r8*cff)))
!>
          ad_fac=ad_Qsea(i)*(1.0_r8+0.378_r8*cff/((PairM-0.378_r8*cff)))
          ad_Qsea(i)=0.0_r8
!>        tl_fac=0.62197_r8*(tl_cff/(PairM-0.378_r8*cff))
!>
          ad_cff=ad_cff+0.62197_r8*ad_fac/(PairM-0.378_r8*cff)
          ad_fac=0.0_r8
!
!  Adjoint Vapor Pressure reduced for salinity
!  (Kraus & Businger, 1994, pp 42).
!
!>        tl_cff=tl_cff*0.98_r8
!>
          ad_cff=ad_cff*0.98_r8
!
!  Adjoint water saturation vapor pressure (mb), using Teten formula.
!
          fac=17.502_r8*TseaC(i)/(240.97_r8+TseaC(i))
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
!>        tl_cff=tl_fac*cff
!>
          ad_fac=ad_cff*cff
          ad_cff=0.0_r8
!>        tl_fac=17.502_r8*tl_TseaC(i)/(240.97_r8+TseaC(i))-            &
!>   &           tl_TseaC(i)*fac/(240.97_r8+TseaC(i))
!>
          ad_TseaC(i)=ad_TseaC(i)+                                      &
     &                17.502_r8*ad_fac/(240.97_r8+TseaC(i))-            &
     &                fac*ad_fac/(240.97_r8+TseaC(i))
          ad_fac=0.0_r8
!
!-----------------------------------------------------------------------
!  Adjoint net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------
!
# ifdef MASKING
!>        tl_LRad(i,j)=tl_LRad(i,j)*rmask(i,j)
!>
          ad_LRad(i,j)=ad_LRad(i,j)*rmask(i,j)
# endif
# if defined LONGWAVE
!
!  Use Berliand (1952) formula to calculate net longwave radiation.
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
!  is assumed constant, but it is a function of latitude varying from
!  1.0 at poles to 0.5 at the Equator).
!
          cff2=TairK(i)*TairK(i)*TairK(i)
!>        tl_LRad(i,j)=-emmiss*StefBo*cff2*4.0_r8*tl_TseaK(i)
!>
          ad_TseaK(i)=ad_TseaK(i)-emmiss*StefBo*cff2*4.0_r8*ad_LRad(i,j)
          ad_LRad(i,j)=0.0_r8

# elif defined LONGWAVE_OUT
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
!>        tl_LRad(i,j)=tl_lrflx(i,j)*Hscale-                            &
!>   &                 4.0_r8*emmiss*StefBo*                            &
!>   &                 tl_TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)
!>
          ad_lrflx(i,j)=ad_lrflx(i,j)+ad_LRad(i,j)*Hscale
          ad_TseaK(i)=ad_TseaK(i)-                                      &
     &                4.0_r8*emmiss*StefBo*ad_LRad(i,j)*                &
     &                TseaK(i)*TseaK(i)*TseaK(i)
          ad_LRad(i,j)=0.0_r8
# else
!>        tl_LRad(i,j)=tl_lrflx(i,j)*Hscale
!>
          ad_lrflx(i,j)=ad_lrflx(i,j)+ad_LRad(i,j)*Hscale
          ad_LRad(i,j)=0.0_r8
# endif
!
!  Initialize.
!
!>        tl_delTc(i)=0.0_r8
!>
          ad_delQc(i)=0.0_r8

!>        tl_delQc(i)=0.0_r8
!>
          ad_delTc(i)=0.0_r8
!>        tl_LHeat(i,j)=tl_lhflx(i,j)*Hscale
!>
          ad_lhflx(i,j)=ad_lhflx(i,j)+ad_LHeat(i,j)*Hscale
          ad_LHeat(i,j)=0.0_r8
!>        tl_SHeat(i,j)=tl_shflx(i,j)*Hscale
!>
          ad_shflx(i,j)=ad_shflx(i,j)+ad_SHeat(i,j)*Hscale
          ad_SHeat(i,j)=0.0_r8
!>        tl_Taux(i,j)=0.0_r8
!>
          ad_Taux(i,j)=0.0_r8
!>        tl_Tauy(i,j)=0.0_r8
!>
          ad_Tauy(i,j)=0.0_r8
!>        tl_Scff(i)=tl_beta(i,j)
!>
          ad_beta(i,j)=ad_beta(i,j)+ad_Scff(i)
          ad_Scff(i)=0.0_r8
!>        tl_Tcff(i)=tl_alpha(i,j)
!>
          ad_alpha(i,j)=ad_alpha(i,j)+ad_Tcff(i)
          ad_Tcff(i)=0.0_r8
!>        tl_rhoSea(i)=tl_rho(i,j,N(ng))
!>
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))+ad_rhoSea(i)
          ad_rhoSea(i)=0.0_r8
!>        tl_TseaK(i)=tl_TseaC(i)
!>
          ad_TseaC(i)=ad_TseaC(i)+ad_TseaK(i)
          ad_TseaK(i)=0.0_r8
!>        tl_TseaC(i)=tl_t(i,j,N(ng),nrhs,itemp)
!>
          ad_t(i,j,N(ng),nrhs,itemp)=ad_t(i,j,N(ng),nrhs,itemp)+        &
     &                               ad_TseaC(i)
          ad_TseaC(i)=0.0_r8
        END DO
      END DO
      RETURN
      END SUBROUTINE ad_bulk_flux_tile

      FUNCTION ad_bulk_psiu (tl_ZoL, ZoL, pi)
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
      real(r8) :: ad_bulk_psiu
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: tl_ZoL, ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Fw, cff, psic, psik, x, y, fac
      real(r8) :: tl_Fw, tl_cff, tl_psic, tl_psik, tl_x, tl_y, tl_fac
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.25_r8
        tl_x=-0.25_r8*15.0_r8*tl_ZoL/((1.0_r8-15.0_r8*ZoL)**0.75_r8)
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))+                             &
     &       LOG(0.5_r8*(1.0_r8+x*x))-                                  &
     &       2.0_r8*ATAN(x)+0.5_r8*pi
        tl_psik=tl_x/(0.5_r8*(1.0_r8+x))+                               &
     &       tl_x*x/(0.5_r8*(1.0_r8+x*x))-                              &
     &       2.0_r8*tl_x/(1.0_r8+x*x)
!     &       2.0_r8*tl_x/SQRT(1.0_r8+x*x)
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-10.15_r8*ZoL)**r3
        tl_y=-r3*10.15_r8*tl_ZoL*(1.0_r8-10.15_r8*ZoL)**(r3-1.0_r8)
        fac=(1.0_r8+2.0_r8*y)/cff
        tl_fac=2.0_r8*tl_y/cff
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN(fac)+pi/cff
        tl_psic=tl_y*(1.0_r8+2.0_r8*y)*1.5_r8/(1.0_r8+y+y*y)-           &
     &       cff*tl_fac/(1.0_r8+fac*fac)
!     &       cff*tl_fac/SQRT(1.0_r8+fac*fac)
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        tl_cff=2.0_r8*tl_ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        tl_Fw=tl_cff/(1.0_r8+cff)-tl_cff*Fw*Fw/cff
        ad_bulk_psiu=(1.0_r8-Fw)*tl_psik-tl_Fw*psik                     &
     &                                         +tl_Fw*psic+Fw*tl_psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        tl_cff=(0.5_r8-SIGN(0.5_r8,0.35_r8*ZoL-50.0))*0.35_r8*tl_ZoL
        fac=EXP(cff)
        tl_fac=tl_cff*fac
        ad_bulk_psiu=-(tl_ZoL+0.6667_r8*tl_ZoL/fac-                     &
     &                 tl_fac*0.6667_r8*(ZoL-14.28_r8)/(fac*fac))
      END IF
      RETURN
      END FUNCTION ad_bulk_psiu

      FUNCTION ad_bulk_psit (tl_ZoL, ZoL, pi)
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
      real(r8) :: ad_bulk_psit
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: tl_ZoL, ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Fw, cff, psic, psik, x, y, fac
      real(r8) :: tl_Fw, tl_cff, tl_psic, tl_psik, tl_x, tl_y, tl_fac
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.5_r8
        IF (x.ne.0.0) THEN
           tl_x=-0.5_r8*15.0_r8*tl_ZoL/x
        ELSE
           tl_x=0.0_r8
        END IF
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))
        tl_psik=tl_x/(0.5_r8*(1.0_r8+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-34.15_r8*ZoL)**r3
        tl_y=-r3*34.15_r8*tl_ZoL*(1.0_r8-34.15_r8*ZoL)**(r3-1.0_r8)
        fac=(1.0_r8+2.0_r8*y)/cff
        tl_fac=2.0_r8*tl_y/cff
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN(fac)+pi/cff
        tl_psic=tl_y*(1.0_r8+2.0_r8*y)*1.5_r8/(1.0_r8+y+y*y)-           &
     &       cff*tl_fac/(1.0_r8+fac*fac)
!     &       cff*tl_fac/SQRT(1.0_r8+fac*fac)
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        tl_cff=2.0_r8*tl_ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        tl_Fw=tl_cff/(1.0_r8+cff)-tl_cff*Fw*Fw/cff
        ad_bulk_psit=(1.0_r8-Fw)*tl_psik-tl_Fw*psik                     &
     &                        +tl_Fw*psic+Fw*tl_psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        tl_cff=(0.5_r8-SIGN(0.5_r8,0.35_r8*ZoL-50.0_r8))*0.35_r8*tl_ZoL
        fac=EXP(cff)
        tl_fac=tl_cff*fac
        ad_bulk_psit=-(3.0_r8*tl_ZoL*(1.0_r8+2.0_r8*ZoL)**0.5_r8+       &
     &            0.6667_r8*tl_ZoL/fac-                                 &
     &               tl_fac*0.6667_r8*(ZoL-14.28_r8)/(fac*fac))
      END IF
      RETURN
      END FUNCTION ad_bulk_psit
#endif
      END MODULE ad_bulk_flux_mod
