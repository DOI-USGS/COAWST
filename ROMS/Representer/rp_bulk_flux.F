#include "cppdefs.h"
      MODULE rp_bulk_flux_mod
#if defined TL_IOMS && defined BULK_FLUXES
!
!svn $Id: rp_bulk_flux.F 795 2016-05-11 01:42:43Z arango $
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
!
      PRIVATE
      PUBLIC  :: rp_bulk_flux
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_bulk_flux (ng, tile)
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
      CALL wclock_on (ng, iRPM, 17)
# endif
      CALL rp_bulk_flux_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng),                                 &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef WET_DRY
     &                        GRID(ng) % rmask_wet,                     &
     &                        GRID(ng) % umask_wet,                     &
     &                        GRID(ng) % vmask_wet,                     &
# endif
     &                        MIXING(ng) % alpha,                       &
     &                        MIXING(ng) % tl_alpha,                    &
     &                        MIXING(ng) % beta,                        &
     &                        MIXING(ng) % tl_beta,                     &
     &                        OCEAN(ng) % rho,                          &
     &                        OCEAN(ng) % tl_rho,                       &
     &                        OCEAN(ng) % t,                            &
     &                        OCEAN(ng) % tl_t,                         &
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
     &                        FORCES(ng) % tl_lhflx,                    &
     &                        FORCES(ng) % lrflx,                       &
     &                        FORCES(ng) % tl_lrflx,                    &
     &                        FORCES(ng) % shflx,                       &
     &                        FORCES(ng) % tl_shflx,                    &
     &                        FORCES(ng) % srflx,                       &
     &                        FORCES(ng) % stflx,                       &
     &                        FORCES(ng) % tl_stflx,                    &
# ifdef EMINUSP
     &                        FORCES(ng) % evap,                        &
     &                        FORCES(ng) % tl_evap,                     &
# endif
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % tl_sustr,                    &
     &                        FORCES(ng) % svstr,                       &
     &                        FORCES(ng) % tl_svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 17)
# endif

      RETURN
      END SUBROUTINE rp_bulk_flux
!
!***********************************************************************
      SUBROUTINE rp_bulk_flux_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs,                               &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef WET_DRY
     &                              rmask_wet, umask_wet, vmask_wet,    &
# endif
     &                              alpha, tl_alpha,                    &
     &                              beta, tl_beta,                      &
     &                              rho, tl_rho, t, tl_t,               &
     &                              Hair, Pair, Tair,                   &
     &                              Uwind, Vwind,                       &
# ifdef CLOUDS
     &                              cloud,                              &
# endif
# ifdef BBL_MODEL_NOT_YET
     &                              Awave, Pwave,                       &
# endif
     &                              rain,                               &
     &                              lhflx, tl_lhflx,                    &
     &                              lrflx, tl_lrflx,                    &
     &                              shflx, tl_shflx,                    &
     &                              srflx,                              &
     &                              stflx, tl_stflx,                    &
# ifdef EMINUSP
     &                              evap, tl_evap,                      &
# endif
     &                              sustr, tl_sustr,                    &
     &                              svstr, tl_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bulk_flux_mod, ONLY : bulk_psiu, bulk_psit
      USE tl_bulk_flux_mod, ONLY : tl_bulk_psiu, tl_bulk_psit
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
      real(r8), intent(in) :: tl_alpha(LBi:,LBj:)
      real(r8), intent(in) :: tl_beta(LBi:,LBj:)
      real(r8), intent(in) :: tl_rho(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: tl_lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: tl_lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: tl_shflx(LBi:,LBj:)
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)

#  ifdef EMINUSP
      real(r8), intent(out) :: evap(LBi:,LBj:)
      real(r8), intent(out) :: tl_evap(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
      real(r8), intent(out) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(out) :: tl_svstr(LBi:,LBj:)
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
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
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
      real(r8), intent(inout) :: tl_lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))

#  ifdef EMINUSP
      real(r8), intent(out) :: evap(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
      integer, parameter :: IterMax = 3

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Cd, Hl, Hlw, Hscale, Hs, Hsr, IER
      real(r8) :: tl_Bf, tl_Cd, tl_Hl, tl_Hlw, tl_Hsr, tl_Hs
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL

      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(r8) :: twopi_inv, wet_bulb
      real(r8) :: tl_wet_bulb, tl_Wspeed, tl_upvel
      real(r8) :: fac, tl_fac, fac1, tl_fac1, fac2, tl_fac2
      real(r8) :: fac3, tl_fac3, tl_cff, ad_upvel
      real(r8) :: adfac, adfac1, adfac2
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
      real(r8) :: tl_Clam,  tl_Hcool, tl_Hsb, tl_Hlb
      real(r8) :: tl_Qbouy, tl_Qcool, tl_lambd
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

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy

      real(r8), dimension(IminS:ImaxS) :: tl_CC
      real(r8), dimension(IminS:ImaxS) :: tl_charn
      real(r8), dimension(IminS:ImaxS) :: tl_Ct
      real(r8), dimension(IminS:ImaxS) :: tl_Cd10
      real(r8), dimension(IminS:ImaxS) :: tl_Ct10
      real(r8), dimension(IminS:ImaxS) :: tl_Cwet
      real(r8), dimension(IminS:ImaxS) :: tl_delQ
      real(r8), dimension(IminS:ImaxS) :: tl_delQc
      real(r8), dimension(IminS:ImaxS) :: tl_delT
      real(r8), dimension(IminS:ImaxS) :: tl_delTc
      real(r8), dimension(IminS:ImaxS) :: tl_delW
      real(r8), dimension(IminS:ImaxS) :: tl_L
      real(r8), dimension(IminS:ImaxS) :: tl_L10
      real(r8), dimension(IminS:ImaxS) :: tl_Q
      real(r8), dimension(IminS:ImaxS) :: tl_Qpsi
      real(r8), dimension(IminS:ImaxS) :: tl_Qsea
      real(r8), dimension(IminS:ImaxS) :: tl_Qstar
      real(r8), dimension(IminS:ImaxS) :: tl_rhoSea
      real(r8), dimension(IminS:ImaxS) :: tl_Ri
      real(r8), dimension(IminS:ImaxS) :: tl_Rr
      real(r8), dimension(IminS:ImaxS) :: tl_Scff
      real(r8), dimension(IminS:ImaxS) :: tl_Tcff
      real(r8), dimension(IminS:ImaxS) :: tl_Tpsi
      real(r8), dimension(IminS:ImaxS) :: tl_TseaC
      real(r8), dimension(IminS:ImaxS) :: tl_TseaK
      real(r8), dimension(IminS:ImaxS) :: tl_Tstar
      real(r8), dimension(IminS:ImaxS) :: tl_u10
      real(r8), dimension(IminS:ImaxS) :: tl_Wgus
      real(r8), dimension(IminS:ImaxS) :: tl_Wpsi
      real(r8), dimension(IminS:ImaxS) :: tl_Wstar
      real(r8), dimension(IminS:ImaxS) :: tl_Zetu
      real(r8), dimension(IminS:ImaxS) :: tl_Zo10
      real(r8), dimension(IminS:ImaxS) :: tl_ZoT10
      real(r8), dimension(IminS:ImaxS) :: tl_ZoL
      real(r8), dimension(IminS:ImaxS) :: tl_ZoQ
      real(r8), dimension(IminS:ImaxS) :: tl_ZoT
      real(r8), dimension(IminS:ImaxS) :: tl_ZoW
      real(r8), dimension(IminS:ImaxS) :: tl_ZWoL

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Tauy

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      twopi_inv=0.5_r8/pi
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          tl_TseaC(i)=tl_t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          tl_TseaK(i)=tl_TseaC(i)+                                      &
#ifdef TL_IOMS
     &                273.16_r8
#endif
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          tl_rhoSea(i)=tl_rho(i,j,N(ng))+                               &
#ifdef TL_IOMS
     &                1000.0_r8
#endif
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
          Tcff(i)=alpha(i,j)
          tl_Tcff(i)=tl_alpha(i,j)
          Scff(i)=beta(i,j)
          tl_Scff(i)=tl_beta(i,j)
!
!  Initialize.
!
          delTc(i)=0.0_r8
          tl_delTc(i)=0.0_r8
          delQc(i)=0.0_r8
          tl_delQc(i)=0.0_r8
          LHeat(i,j)=lhflx(i,j)*Hscale
          tl_LHeat(i,j)=tl_lhflx(i,j)*Hscale
          SHeat(i,j)=shflx(i,j)*Hscale
          tl_SHeat(i,j)=tl_shflx(i,j)*Hscale
          Taur=0.0_r8
          Taux(i,j)=0.0_r8
          tl_Taux(i,j)=0.0_r8
          Tauy(i,j)=0.0_r8
          tl_Tauy(i,j)=0.0_r8
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
          tl_LRad(i,j)=-emmiss*StefBo*cff2*4.0_r8*tl_TseaK(i)-          &
#ifdef TL_IOMS
     &                  emmiss*StefBo*                                  &
     &                  (cff1*(0.39_r8-0.05_r8*SQRT(vap_p))*            &
     &                        (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))- &
     &                   cff2*4.0_r8*TairK(i))
#endif
# elif defined LONGWAVE_OUT
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          LRad(i,j)=lrflx(i,j)*Hscale-                                  &
     &              emmiss*StefBo*TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)
          tl_LRad(i,j)=tl_lrflx(i,j)*Hscale                             &
     &                 -4.0_r8*emmiss*StefBo*tl_TseaK(i)*               &
     &                             TseaK(i)*TseaK(i)*TseaK(i)+          &
#ifdef TL_IOMS
     &         3.0_r8*emmiss*StefBo*TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)
#endif

# else
          LRad(i,j)=lrflx(i,j)*Hscale
#  ifdef TL_IOMS
          tl_LRad(i,j)=lrflx(i,j)*Hscale
#  else
          tl_LRad(i,j)=tl_lrflx(i,j)*Hscale
#  endif
# endif
# ifdef MASKING
          LRad(i,j)=LRad(i,j)*rmask(i,j)
          tl_LRad(i,j)=tl_LRad(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          LRad(i,j)=LRad(i,j)*rmask_wet(i,j)
          tl_LRad(i,j)=tl_LRad(i,j)*rmask_wet(i,j)
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
          fac=17.502_r8*TseaC(i)/(240.97_r8+TseaC(i))
          tl_fac=17.502_r8*tl_TseaC(i)/(240.97_r8+TseaC(i))-            &
     &           tl_TseaC(i)*fac/(240.97_r8+TseaC(i))+                  &
#ifdef TL_IOMS
     &           fac*TseaC(i)/(240.97_r8+TseaC(i))
#endif

!>        cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*EXP(fac)
!>
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
          tl_cff=tl_fac*cff+                                            &
#ifdef TL_IOMS
     &           (1.0_r8-fac)*cff
#endif
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
          cff=cff*0.98_r8
          tl_cff=tl_cff*0.98_r8
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
          Qsea(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
          tl_Qsea(i)=0.62197_r8*(tl_cff/(PairM-0.378_r8*cff))           &
     &               *(1.0_r8+0.378_r8*cff/((PairM-0.378_r8*cff)))-     &
#ifdef TL_IOMS
     &                       0.378_r8*Qsea(i)*cff/(PairM-0.378_r8*cff)
#endif
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
          tl_Hlv(i,j)=-0.00237_r8*tl_TseaC(i)*1.0E+6_r8+                &
#ifdef TL_IOMS
     &                 2.501E+6_r8
#endif
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus(i)=0.5_r8
#ifdef TL_IOMS
          tl_Wgus(i)=0.5_r8
#else
          tl_Wgus(i)=0.0_r8
#endif
!
!  Note Wgus is a constant at this point.
!
          delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
#ifdef TL_IOMS
          tl_delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
#else
          tl_delW(i)=0.0_r8
#endif
          delQ(i)=Qsea(i)-Q(i)
          tl_delQ(i)=tl_Qsea(i)-                                        &
#ifdef TL_IOMS
     &                        Q(i)
#endif
          delT(i)=TseaC(i)-TairC(i)
          tl_delT(i)=tl_TseaC(i)-                                       &
#ifdef TL_IOMS
     &               TairC(i)
#endif
!
!  Neutral coefficients.
!
          ZoW(i)=0.0001_r8
#ifdef TL_IOMS
          tl_ZoW(i)=0.0001_r8
#else
          tl_ZoW(i)=0.0_r8
#endif
          u10(i)=delW(i)*LOG(10.0_r8/ZoW(i))/LOG(blk_ZW(ng)/ZoW(i))
!
!  Note: ZoW is a constant at this point
!
          tl_u10(i)=tl_delW(i)*LOG(10.0_r8/ZoW(i))/                     &
     &                         LOG(blk_ZW(ng)/ZoW(i))
          Wstar(i)=0.035_r8*u10(i)
          tl_Wstar(i)=0.035_r8*tl_u10(i)
          Zo10(i)=0.011_r8*Wstar(i)*Wstar(i)/g+                         &
     &            0.11_r8*VisAir(i)/Wstar(i)
          tl_Zo10(i)=0.022_r8*tl_Wstar(i)*Wstar(i)/g-                   &
     &            tl_Wstar(i)*0.11_r8*VisAir(i)/(Wstar(i)*Wstar(i))-    &
#ifdef TL_IOMS
     &          0.011_r8*Wstar(i)*Wstar(i)/g+0.22_r8*VisAir(i)/Wstar(i)
#endif
          fac=LOG(10.0_r8/Zo10(i))
          tl_fac=-tl_Zo10(i)/Zo10(i)+                                   &
#ifdef TL_IOMS
     &           fac+1.0_r8
#endif
!>          Cd10(i)=(vonKar/fac)**2
          Cd10(i)=(vonKar/LOG(10.0_r8/Zo10(i)))**2
          tl_Cd10(i)=-2.0_r8*tl_fac*Cd10(i)/fac+                        &
#ifdef TL_IOMS
     &               3.0_r8*Cd10(i)
#endif
          Ch10(i)=0.00115_r8
          Ct10(i)=Ch10(i)/sqrt(Cd10(i))
          tl_Ct10(i)=-0.5_r8*tl_Cd10(i)*Ct10(i)/Cd10(i)+                &
#ifdef TL_IOMS
     &               1.5_r8*Ct10(i)
#endif
          fac=vonKar/Ct10(i)
          tl_fac=-tl_Ct10(i)*fac/Ct10(i)+                               &
#ifdef TL_IOMS
     &               2.0_r8*fac
#endif
!>          ZoT10(i)=10.0_r8/EXP(fac)
          ZoT10(i)=10.0_r8/EXP(vonKar/Ct10(i))
          tl_ZoT10(i)=-tl_fac*ZoT10(i)+                                 &
#ifdef TL_IOMS
     &               (1.0_r8+fac)*ZoT10(i)
#endif
          fac=LOG(blk_ZW(ng)/Zo10(i))
          tl_fac=-tl_Zo10(i)/Zo10(i)+                                   &
#ifdef TL_IOMS
     &             fac+1.0_r8
#endif
!>          Cd=(vonKar/fac)**2
          Cd=(vonKar/LOG(blk_ZW(ng)/Zo10(i)))**2
          tl_Cd=-2.0_r8*tl_fac*Cd/fac+                                  &
#ifdef TL_IOMS
     &               3.0_r8*Cd
#endif
!
!  Compute Richardson number.
!
          fac=LOG(blk_ZT(ng)/ZoT10(i))
          tl_fac=-tl_ZoT10(i)/ZoT10(i)+                                 &
#ifdef TL_IOMS
     &             fac+1.0_r8
#endif
!>          Ct(i)=vonKar/fac                        ! T transfer coefficient
          Ct(i)=vonKar/LOG(blk_ZT(ng)/ZoT10(i))     ! T transfer coefficient
          tl_Ct(i)=-tl_fac*Ct(i)/fac+                                   &
#ifdef TL_IOMS
     &               2.0_r8*Ct(i)
#endif
          CC(i)=vonKar*Ct(i)/Cd
          tl_CC(i)=vonKar*tl_Ct(i)/Cd-tl_Cd*CC(i)/Cd+                   &
#ifdef TL_IOMS
     &               CC(i)
#endif
          delTc(i)=0.0_r8
          tl_delTc(i)=0.0_r8
# ifdef COOL_SKIN
          delTc(i)=blk_dter
#   ifdef TL_IOMS
          tl_delTc(i)=blk_dter
#   else
          tl_delTc(i)=0.0_r8
#   endif
# endif
          Ribcu(i)=-blk_ZW(ng)/(blk_Zabl*0.004_r8*blk_beta**3)
          fac=1.0/(TairK(i)*delW(i)*delW(i))
!
!  Note delW is a constant at this point so we do not need tl_fac.
!
!>        Ri(i)=-fac*g*blk_ZW(ng)*((delT(i)-delTc(i))+                &
!>   &                              0.61_r8*TairK(i)*delQ(i))
!>
          Ri(i)=-g*blk_ZW(ng)*((delT(i)-delTc(i))+                      &
     &                          0.61_r8*TairK(i)*delQ(i))/              &
     &          (TairK(i)*delW(i)*delW(i))
          tl_Ri(i)=-fac*g*blk_ZW(ng)*((tl_delT(i)-tl_delTc(i))+         &
     &                                 0.61_r8*TairK(i)*tl_delQ(i))
          IF (Ri(i).lt.0.0_r8) THEN
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+Ri(i)/Ribcu(i))       ! Unstable
            tl_Zetu(i)=(tl_CC(i)*Ri(i)+CC(i)*tl_Ri(i))/                 &
     &                                 (1.0_r8+Ri(i)/Ribcu(i))-         &
     &              (tl_Ri(i)/Ribcu(i))*Zetu(i)/(1.0_r8+Ri(i)/Ribcu(i))-&
#ifdef TL_IOMS
                      Zetu(i)/(1.0_r8+Ri(i)/Ribcu(i))
#endif
          ELSE
            fac=3.0_r8*Ri(i)/CC(i)
            tl_fac=3.0_r8*tl_Ri(i)/CC(i)-tl_CC(i)*fac/CC(i)+            &
#ifdef TL_IOMS
                      fac
#endif
!>          Zetu(i)=CC(i)*Ri(i)/(1.0_r8+fac)                  ! Stable
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+3.0_r8*Ri(i)/CC(i))   ! Stable
            tl_Zetu(i)=(tl_CC(i)*Ri(i)+CC(i)*tl_Ri(i))/(1.0_r8+fac)     &
     &                -tl_fac*Zetu(i)/(1.0_r8+fac)-                     &
#ifdef TL_IOMS
                      Zetu(i)/(1.0_r8+fac)
#endif
          END IF
          L10(i)=blk_ZW(ng)/Zetu(i)
          tl_L10(i)=-L10(i)*L10(i)*tl_Zetu(i)/blk_ZW(ng)+               &
#ifdef TL_IOMS
     &                 2.0_r8*L10(i)
#endif
!
!  First guesses for Monon-Obukhov similarity scales.
!
          fac1=blk_ZW(ng)/L10(i)
          tl_fac1=-tl_L10(i)*fac1/L10(i)+                               &
#ifdef TL_IOMS
     &            2.0_r8*fac1
#endif
          fac2=LOG(blk_ZW(ng)/Zo10(i))
          tl_fac2=-tl_Zo10(i)/Zo10(i)+                                  &
#ifdef TL_IOMS
     &            fac2+1.0_r8
#endif
          fac3=bulk_psiu(fac1,pi)
#ifdef TL_IOMS
!          tl_fac3=rp_bulk_psiu(tl_fac1,fac1,pi)
          tl_fac3=bulk_psiu(fac1,pi)+                                   &
     &            tl_bulk_psiu(tl_fac1,fac1,pi)-                        &
     &            tl_bulk_psiu(fac1,fac1,pi)
#else
          tl_fac3=tl_bulk_psiu(tl_fac1,fac1,pi)
#endif
!>          Wstar(i)=delW(i)*vonKar/(fac2-bulk_psiu(fac1,pi))
!>          Wstar(i)=delW(i)*vonKar/(fac2-fac3))
          Wstar(i)=delW(i)*vonKar/(LOG(blk_ZW(ng)/Zo10(i))-             &
     &                             bulk_psiu(blk_ZW(ng)/L10(i),pi))
          tl_Wstar(i)=-(tl_fac2-tl_fac3)*Wstar(i)                       &
     &                                 /(fac2-fac3)+                    &
#ifdef TL_IOMS
     &                 2.0_r8*Wstar(i)
#endif
          fac1=blk_ZT(ng)/L10(i)
          tl_fac1=-tl_L10(i)*fac1/L10(i)+                               &
#ifdef TL_IOMS
     &            2.0_r8*fac1
#endif
          fac2=LOG(blk_ZT(ng)/ZoT10(i))
          tl_fac2=-tl_ZoT10(i)/ZoT10(i)+                                &
#ifdef TL_IOMS
     &            fac2+1.0_r8
#endif
          fac3=bulk_psit(fac1,pi)
#ifdef TL_IOMS
!          tl_fac3=rp_bulk_psit(tl_fac1,fac1,pi)
          tl_fac3=bulk_psit(fac1,pi)+                                   &
     &            tl_bulk_psit(tl_fac1,fac1,pi)-                        &
     &            tl_bulk_psit(fac1,fac1,pi)
#else
          tl_fac3=tl_bulk_psit(tl_fac1,fac1,pi)
#endif
!>          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
!>     &             (fac2-bulk_psit(fac1,pi))
!>          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
!>     &             (fac2-fac3)
          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                          &
     &             (LOG(blk_ZT(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZT(ng)/L10(i),pi))
          tl_Tstar(i)=-(tl_delT(i)-tl_delTc(i))*vonKar/                 &
     &             (fac2-fac3)-                                         &
     &             (tl_fac2-tl_fac3)*Tstar(i)/                          &
     &             (fac2-fac3)+                                         &
#ifdef TL_IOMS
     &            Tstar(i)
#endif
          fac1=blk_ZQ(ng)/L10(i)
          tl_fac1=-tl_L10(i)*fac1/L10(i)+                               &
#ifdef TL_IOMS
     &            2.0_r8*fac1
#endif
          fac2=LOG(blk_ZQ(ng)/ZoT10(i))
          tl_fac2=-tl_ZoT10(i)/ZoT10(i)+                                &
#ifdef TL_IOMS
     &            fac2+1.0_r8
#endif
          fac3=bulk_psit(fac1,pi)
#ifdef TL_IOMS
!          tl_fac3=rp_bulk_psit(tl_fac1,fac1,pi)
          tl_fac3=bulk_psit(fac1,pi)+                                   &
     &            tl_bulk_psit(tl_fac1,fac1,pi)-                        &
     &            tl_bulk_psit(fac1,fac1,pi)
#else
          tl_fac3=tl_bulk_psit(tl_fac1,fac1,pi)
#endif
!>          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
!>     &             (fac2-bulk_psit(fac1,pi))
!>          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
!>     &             (fac2-fac3)
          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                          &
     &             (LOG(blk_ZQ(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZQ(ng)/L10(i),pi))
          tl_Qstar(i)=-(tl_delQ(i)-tl_delQc(i))*vonKar/                 &
     &             (fac2-fac3)-                                         &
     &             (tl_fac2-tl_fac3)*Qstar(i)/                          &
     &             (fac2-fac3)+                                         &
#ifdef TL_IOMS
     &            Qstar(i)
#endif
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          IF (delW(i).gt.18.0_r8) THEN
            charn(i)=0.018_r8
#ifdef TL_IOMS
            tl_charn(i)=0.018_r8
#else
            tl_charn(i)=0.0_r8
#endif
          ELSE IF ((10.0_r8.lt.delW(i)).and.(delW(i).le.18.0_r8)) THEN
            charn(i)=0.011_r8+                                          &
     &               0.125_r8*(0.018_r8-0.011_r8)*(delW(i)-10._r8)
#ifdef TL_IOMS
            tl_charn(i)=0.011_r8+                                       &
     &               0.125_r8*(0.018_r8-0.011_r8)*(delW(i)-10._r8)
#else
            tl_charn(i)=0.0_r8
#endif
          ELSE
            charn(i)=0.011_r8
#ifdef TL_IOMS
            tl_charn(i)=0.011_r8
#else
            tl_charn(i)=0.0_r8
#endif
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
            tl_ZoW(i)=(25._r8/pi)*Lwave(i)*4.5_r8*tl_Wstar(i)*          &
     &              (Wstar(i)/Cwave(i))**3.5/Cwave(i)-                  &
     &             tl_Wstar(i)*0.11_r8*VisAir(i)/                       &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))-    &
#    ifdef TL_IOMS
     &             3.5_r8*(25._r8/pi)*Lwave(i)*(Wstar(i)/Cwave(i))**4.5+&
     &             0.11_r8*VisAir(i)*(2.0_r8*Wstar(i)+eps)/             &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
#    endif
#  else
            ZoW(i)=1200._r8*Awave(i,j)*(Awave(i,j)/Lwave(i))**4.5+      &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
            tl_ZoW(i)=-tl_Wstar(i)*0.11_r8*VisAir(i)/                   &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))+    &
#    ifdef TL_IOMS
     &              0.11_r8*VisAir(i)*(2.0_r8*Wstar(i)+eps)/            &
     &                              ((Wstar(i)+eps)*(Wstar(i)+eps))
#    endif
#  endif
# else
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
            tl_ZoW(i)=2.0_r8*charn(i)*tl_Wstar(i)*Wstar(i)/g-           &
     &             tl_Wstar(i)*0.11_r8*VisAir(i)/                       &
     &                           ((Wstar(i)+eps)*(Wstar(i)+eps))-       &
#  ifdef TL_IOMS
     &                charn(i)*Wstar(i)*Wstar(i)/g+                     &
     &               0.11_r8*VisAir(i)*(2.0_r8*Wstar(i)+eps)/           &
     &                           ((Wstar(i)+eps)*(Wstar(i)+eps))
#  endif
# endif
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
            tl_Rr(i)=(tl_ZoW(i)*Wstar(i)+ZoW(i)*tl_Wstar(i))/VisAir(i)- &
# ifdef TL_IOMS
     &            Rr(i)
# endif
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6)
            tl_ZoQ(i)=                                                  &
     &         -(0.5_r8-SIGN(0.5_r8,5.5e-5_r8/Rr(i)**0.6-1.15e-4_r8))   &
     &          *0.6_r8*5.5e-5_r8*tl_Rr(i)/Rr(i)**1.6+                  &
# ifdef TL_IOMS
     &          (0.5_r8-SIGN(0.5_r8,5.5e-5_r8/Rr(i)**0.6-1.15e-4_r8))   &
     &          *5.5e-5_r8*1.6_r8/Rr(i)**0.6+                           &
     &          (0.5_r8+SIGN(0.5_r8,5.5e-5_r8/Rr(i)**0.6-1.15e-4_r8))*  &
     &          1.15e-4_r8
# endif
            ZoT(i)=ZoQ(i)
            tl_ZoT(i)=tl_ZoQ(i)
            ZoL(i)=vonKar*g*blk_ZW(ng)*                                 &
     &             (Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                     &
     &                        0.61_r8*TairK(i)*Qstar(i))/               &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
            tl_ZoL(i)=vonKar*g*blk_ZW(ng)*                              &
     &             (tl_Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                  &
     &                        0.61_r8*TairK(i)*tl_Qstar(i))/            &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)-                          &
     &             2.0_r8*TairK(i)*tl_Wstar(i)*Wstar(i)*                &
     &             (1.0_r8+0.61_r8*Q(i))*ZoL(i)/                        &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)+                          &
# ifdef TL_IOMS
     &           2.0_r8*ZoL(i)*TairK(i)*Wstar(i)*Wstar(i)*              &
     &                                           (1.0_r8+0.61_r8*Q(i))/ &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
# endif

            L(i)=blk_ZW(ng)/(ZoL(i)+eps)
            tl_L(i)=-tl_ZoL(i)*L(i)/(ZoL(i)+eps)+                       &
# ifdef TL_IOMS
     &           L(i)*(2.0_r8*ZoL(i)+eps)/(ZoL(i)+eps)
# endif
!
!  Evaluate stability functions at Z/L.
!
            Wpsi(i)=bulk_psiu(ZoL(i),pi)
!            tl_Wpsi(i)=rp_bulk_psiu(tl_ZoL(i),ZoL(i),pi)
# ifdef TL_IOMS
            tl_Wpsi(i)=bulk_psiu(ZoL(i),pi)+                            &
     &                 tl_bulk_psiu(tl_ZoL(i),ZoL(i),pi)-               &
     &                 tl_bulk_psiu(ZoL(i),ZoL(i),pi)
# else
            tl_Wpsi(i)=tl_bulk_psiu(tl_ZoL(i),ZoL(i),pi)
# endif
            fac=blk_ZT(ng)/L(i)
            tl_fac=-tl_L(i)*fac/L(i)+                                   &
# ifdef TL_IOMS
     &              2.0_r8*fac
# endif
!>            Tpsi(i)=bulk_psit(fac,pi)
            Tpsi(i)=bulk_psit(blk_ZT(ng)/L(i),pi)
!            tl_Tpsi(i)=rp_bulk_psit(tl_fac,fac,pi)
# ifdef TL_IOMS
            tl_Tpsi(i)=bulk_psit(fac,pi)+                               &
     &                 tl_bulk_psit(tl_fac,fac,pi)-                     &
     &                 tl_bulk_psit(fac,fac,pi)
# else
            tl_Tpsi(i)=tl_bulk_psit(tl_fac,fac,pi)
# endif
            fac=blk_ZQ(ng)/L(i)
            tl_fac=-tl_L(i)*fac/L(i)+                                   &
# ifdef TL_IOMS
     &              2.0_r8*fac
# endif
!>            Qpsi(i)=bulk_psit(fac,pi)
            Qpsi(i)=bulk_psit(blk_ZQ(ng)/L(i),pi)
!            tl_Qpsi(i)=rp_bulk_psit(tl_fac,fac,pi)
# ifdef TL_IOMS
            tl_Qpsi(i)=bulk_psit(fac,pi)+                               &
     &                 tl_bulk_psit(tl_fac,fac,pi)-                     &
     &                 tl_bulk_psit(fac,fac,pi)
# else
            tl_Qpsi(i)=tl_bulk_psit(tl_fac,fac,pi)
# endif
# ifdef COOL_SKIN
            Cwet(i)=0.622_r8*Hlv(i,j)*Qsea(i)/                          &
     &              (blk_Rgas*TseaK(i)*TseaK(i))
            tl_Cwet(i)=0.622_r8*(tl_Hlv(i,j)*Qsea(i)+                   &
     &                                           Hlv(i,j)*tl_Qsea(i))/  &
     &              (blk_Rgas*TseaK(i)*TseaK(i))-                       &
     &               2.0_r8*blk_Rgas*tl_TseaK(i)*TseaK(i)*Cwet(i)/      &
     &              (blk_Rgas*TseaK(i)*TseaK(i))+                       &
# ifdef TL_IOMS
     &              Cwet(i)
# endif
            delQc(i)=Cwet(i)*delTc(i)
            tl_delQc(i)=tl_Cwet(i)*delTc(i)+Cwet(i)*tl_delTc(i)-        &
#  ifdef TL_IOMS
     &               delQc(i)
#  endif
# endif
!
!  Compute wind scaling parameters, Wstar.
!
            fac1=blk_ZW(ng)/ZoW(i)
            tl_fac1=-tl_ZoW(i)*fac1/ZoW(i)+                             &
#  ifdef TL_IOMS
     &               2.0_r8*fac1
#  endif
            fac2=LOG(fac1)
            tl_fac2=tl_fac1/fac1+                                       &
#  ifdef TL_IOMS
     &               fac2-1.0_r8
#  endif
            fac3=delW(i)*vonKar/(fac2-Wpsi(i))
            tl_fac3=tl_delW(i)*vonKar/(fac2-Wpsi(i))-                   &
     &             (tl_fac2-tl_Wpsi(i))*fac3/(fac2-Wpsi(i))+            &
#  ifdef TL_IOMS
     &               fac3
#  endif
!>            Wstar(i)=MAX(eps,fac3)
            Wstar(i)=MAX(eps,delW(i)*vonKar/                            &
     &               (LOG(blk_ZW(ng)/ZoW(i))-Wpsi(i)))

            tl_Wstar(i)=(0.5_r8-SIGN(0.5_r8,eps-fac3))*tl_fac3+         &
#  ifdef TL_IOMS
     &                  (0.5_r8+SIGN(0.5_r8,eps-fac3))*eps
#  endif
            fac1=blk_ZT(ng)/ZoT(i)
            tl_fac1=-tl_ZoT(i)*fac1/ZoT(i)+                             &
#  ifdef TL_IOMS
     &               2.0_r8*fac1
#  endif
            fac2=LOG(fac1)
            tl_fac2=tl_fac1/fac1+                                       &
#  ifdef TL_IOMS
     &               fac2-1.0_r8
#  endif
!>            Tstar(i)=-(delT(i)-delTc(i))*vonKar/(fac2-Tpsi(i))
            Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
     &               (LOG(blk_ZT(ng)/ZoT(i))-Tpsi(i))
            tl_Tstar(i)=-(tl_delT(i)-tl_delTc(i))*vonKar/(fac2-Tpsi(i))-&
     &           (tl_fac2-tl_Tpsi(i))*Tstar(i)/(fac2-Tpsi(i))+          &
#  ifdef TL_IOMS
     &                   Tstar(i)
#  endif
            fac1=blk_ZQ(ng)/ZoQ(i)
            tl_fac1=-tl_ZoQ(i)*fac1/ZoQ(i)+                             &
#  ifdef TL_IOMS
     &               2.0_r8*fac1
#  endif
            fac2=LOG(fac1)
            tl_fac2=tl_fac1/fac1+                                       &
#  ifdef TL_IOMS
     &               fac2-1.0_r8
#  endif
!>            Qstar(i)=-(delQ(i)-delQc(i))*vonKar/(fac2-Qpsi(i))
            Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
     &               (LOG(blk_ZQ(ng)/ZoQ(i))-Qpsi(i))
            tl_Qstar(i)=-(tl_delQ(i)-tl_delQc(i))*vonKar/(fac2-Qpsi(i))-&
     &                    (tl_fac2-tl_Qpsi(i))*Qstar(i)/(fac2-Qpsi(i))+ &
#  ifdef TL_IOMS
     &                   Qstar(i)
#  endif
!
!  Compute gustiness in wind speed.
!
            Bf=-g/TairK(i)*                                             &
     &         Wstar(i)*(Tstar(i)+0.61_r8*TairK(i)*Qstar(i))
            tl_Bf=-g/TairK(i)*(                                         &
     &         tl_Wstar(i)*(Tstar(i)+0.61_r8*TairK(i)*Qstar(i))+        &
     &         Wstar(i)*(tl_Tstar(i)+0.61_r8*TairK(i)*tl_Qstar(i)) )-   &
#ifdef TL_IOMS
     &                  Bf
#endif
            IF (Bf.gt.0.0_r8) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
              tl_Wgus(i)=r3*blk_beta*tl_Bf*blk_Zabl*                    &
     &                                      (Bf*blk_Zabl)**(r3-1.0_r8)+ &
#ifdef TL_IOMS
     &                    (1.0_r8-r3)*Wgus(i)
#endif
            ELSE
              Wgus(i)=0.2_r8
#ifdef TL_IOMS
              tl_Wgus(i)=0.2_r8
#else
              tl_Wgus(i)=0.0_r8
#endif
            END IF
            delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
            IF (delW(i).ne.0.0_r8)THEN
               tl_delW(i)=tl_Wgus(i)*Wgus(i)/delW(i)+                   &
#ifdef TL_IOMS
     &                    Wmag(i)*Wmag(i)/delW(i)
#endif
            ELSE
               tl_delW(i)=0.0_r8
            END IF
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
            tl_Clam=48.0_r8*g*blk_Cpw*tl_rhoSea(i)*blk_visw*            &
     &                             (rhoSea(i)*blk_visw)**2/             &
     &           (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))-                 &
#ifdef TL_IOMS
     &                    2.0_r8*Clam
#endif
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
            Hcool=0.001_r8
!
!  Backgound sensible and latent heat.
!
            Hsb=-rhoAir(i)*blk_Cpa*Wstar(i)*Tstar(i)
            tl_Hsb=-rhoAir(i)*blk_Cpa*(tl_Wstar(i)*Tstar(i)+            &
     &                                 Wstar(i)*tl_Tstar(i))-           &
#ifdef TL_IOMS
     &                    Hsb
#endif
            Hlb=-rhoAir(i)*Hlv(i,j)*Wstar(i)*Qstar(i)
            tl_Hlb=-rhoAir(i)*(tl_Hlv(i,j)*Wstar(i)*Qstar(i)+           &
     &                         Hlv(i,j)*(tl_Wstar(i)*Qstar(i)+          &
     &                                 Wstar(i)*tl_Qstar(i)))-          &
#ifdef TL_IOMS
     &                    2.0_r8*Hlb
#endif
!
!  Mean absoption in cool-skin layer.
!
            Fc=0.065_r8+11.0_r8*Hcool-                                  &
     &         (1.0_r8-EXP(-Hcool*1250.0_r8))*6.6E-5_r8/Hcool
!
!  Total cooling at the interface.
!
            Qcool=LRad(i,j)+Hsb+Hlb-SRad(i,j)*Fc
            tl_Qcool=tl_LRad(i,j)+tl_Hsb+tl_Hlb-                        &
#ifdef TL_IOMS
     &                     SRad(i,j)*Fc
#endif
            Qbouy=Tcff(i)*Qcool+Scff(i)*Hlb*blk_Cpw/Hlv(i,j)
            tl_Qbouy=tl_Tcff(i)*Qcool+Tcff(i)*tl_Qcool+                 &
     &                  (tl_Scff(i)*Hlb*blk_Cpw+                        &
     &                   Scff(i)*tl_Hlb*blk_Cpw)/Hlv(i,j)-              &
     &           tl_Hlv(i,j)*Scff(i)*Hlb*blk_Cpw/(Hlv(i,j)*Hlv(i,j))-   &
#ifdef TL_IOMS
     &                    Tcff(i)*Qcool
#endif
!
!  Compute temperature and moisture change.
!
            IF ((Qcool.gt.0.0_r8).and.(Qbouy.gt.0.0_r8)) THEN
              fac1=(Wstar(i)+eps)**4
              tl_fac1=4.0_r8*tl_Wstar(i)*(Wstar(i)+eps)**3+             &
#ifdef TL_IOMS
     &                (eps-3.0_r8*Wstar(i))*(Wstar(i)+eps)**3
#endif
              fac2=Clam*Qbouy
              tl_fac2=tl_Clam*Qbouy+Clam*tl_Qbouy-                      &
#ifdef TL_IOMS
     &                fac2
#endif
              fac3=(fac2/fac1)**0.75_r8
              tl_fac3=0.75_r8*(fac2/fac1)**(-0.25_r8)*                  &
     &              (tl_fac2/fac1-tl_fac1*fac2/(fac1*fac1))+            &
#ifdef TL_IOMS
     &                fac3
#endif
!>              lambd=6.0_r8/(1.0_r8+fac3)**r3
            lambd=6.0_r8/                                               &
     &              (1.0_r8+(Clam*Qbouy/(Wstar(i)+eps)**4)**0.75_r8)**r3
              tl_lambd=-r3*6.0_r8*tl_fac3/(1.0_r8+fac3)**(r3+1.0_r8)+   &
#ifdef TL_IOMS
     &              lambd*(1.0_r8+(1.0_r8+r3)*fac3)/(1.0_r8+fac3)
#endif
              fac1=SQRT(rhoAir(i)/rhoSea(i))
              IF (fac1.ne.0.0_r8) THEN
                  tl_fac1=-0.5_r8*tl_rhoSea(i)*fac1/rhoSea(i)+          &
#ifdef TL_IOMS
     &                   1.5_r8*fac1
#endif
              ELSE
                  tl_fac1=0.0_r8
              END IF
              fac2=fac1*Wstar(i)+eps
              tl_fac2=tl_fac1*Wstar(i)+fac1*tl_Wstar(i)-                &
#ifdef TL_IOMS
     &                   fac1*Wstar(i)+eps
#endif
!>              Hcool=lambd*blk_visw/fac2
            Hcool=lambd*blk_visw/(SQRT(rhoAir(i)/rhoSea(i))*            &
     &            Wstar(i)+eps)
              tl_Hcool=tl_lambd*blk_visw/fac2-tl_fac2*Hcool/fac2+       &
#ifdef TL_IOMS
     &                   Hcool
#endif
              delTc(i)=Qcool*Hcool/blk_tcw
              tl_delTc(i)=(tl_Qcool*Hcool+Qcool*tl_Hcool)/blk_tcw-      &
#ifdef TL_IOMS
     &                   delTc(i)
#endif
            ELSE
              delTc(i)=0.0_r8
              tl_delTc(i)=0.0_r8
            END IF
            delQc(i)=Cwet(i)*delTc(i)
            tl_delQc(i)=tl_Cwet(i)*delTc(i)+Cwet(i)*tl_delTc(i)-        &
#ifdef TL_IOMS
     &                   delQc(i)
#endif
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
          IF (Wspeed.ne.0.0_r8) THEN
            tl_Wspeed=tl_Wgus(i)*Wgus(i)/Wspeed+                        &
#ifdef TL_IOMS
     &                   Wmag(i)*Wmag(i)/Wspeed
#endif
          ELSE
            tl_Wspeed=0.0_r8
          END IF
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
          tl_Cd=2.0_r8*(tl_Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)        &
     &          -tl_Wspeed*Wspeed*Cd/(Wspeed*Wspeed+eps))+              &
#ifdef TL_IOMS
     &           Cd*(Wspeed*Wspeed-eps)/(Wspeed*Wspeed+eps)
#endif
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
            Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
            tl_Hs=-blk_Cpa*rhoAir(i)*(tl_Wstar(i)*Tstar(i)+             &
     &                                Wstar(i)*tl_Tstar(i))-            &
#ifdef TL_IOMS
     &                        Hs
#endif
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
            diffw=2.11E-5_r8*(TairK(i)/273.16_r8)**1.94_r8
            diffh=0.02411_r8*(1.0_r8+TairC(i)*                          &
     &                        (3.309E-3_r8-1.44E-6_r8*TairC(i)))/       &
     &            (rhoAir(i)*blk_Cpa)
            cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
            tl_cff=Qair(i)*tl_Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
            fac=0.622_r8*(cff*Hlv(i,j)*diffw)/(blk_Cpa*diffh)
            tl_fac=0.622_r8*diffw*(tl_cff*Hlv(i,j)+cff*tl_Hlv(i,j))/    &
     &                                      (blk_Cpa*diffh)-            &
#ifdef TL_IOMS
     &                        fac
#endif
!>            wet_bulb=1.0_r8/(1.0_r8+fac)
            wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/      &
     &                                       (blk_Cpa*diffh))
            tl_wet_bulb=-tl_fac*wet_bulb*wet_bulb+                      &
#ifdef TL_IOMS
     &                  wet_bulb*wet_bulb*(1.0_r8+2.0_r8*fac)
#endif
            Hsr=rain(i,j)*wet_bulb*blk_Cpw*                             &
     &          ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
            tl_Hsr=Hsr*tl_wet_bulb/wet_bulb+                            &
     &          rain(i,j)*wet_bulb*blk_Cpw*                             &
     &          (tl_TseaC(i)+tl_Qsea(i)*Hlv(i,j)/blk_Cpa+               &
     &          (Qsea(i)-Q(i))*tl_Hlv(i,j)/blk_Cpa)-                    &
#ifdef TL_IOMS
     &         rain(i,j)*wet_bulb*blk_Cpw*                              &
     &            (TseaC(i)+(2.0_r8*Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
#endif
            SHeat(i,j)=(Hs+Hsr)
            tl_SHeat(i,j)=(tl_Hs+tl_Hsr)
# ifdef MASKING
            SHeat(i,j)=SHeat(i,j)*rmask(i,j)
            tl_SHeat(i,j)=tl_SHeat(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
            SHeat(i,j)=SHeat(i,j)*rmask_wet(i,j)
            tl_SHeat(i,j)=tl_SHeat(i,j)*rmask_wet(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
            Hl=-Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)
            tl_Hl=-tl_Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)-             &
     &             Hlv(i,j)*rhoAir(i)*(tl_Wstar(i)*Qstar(i)+            &
     &                                     Wstar(i)*tl_Qstar(i) )-      &
#ifdef TL_IOMS
     &                      2.0_r8*Hl
#endif
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
            upvel=-1.61_r8*Wstar(i)*Qstar(i)-                           &
     &            (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
            tl_upvel=-1.61_r8*                                          &
     &            (tl_Wstar(i)*Qstar(i)+Wstar(i)*tl_Qstar(i))-          &
     &            (1.0_r8+1.61_r8*Q(i))*(tl_Wstar(i)*Tstar(i)+          &
     &                                Wstar(i)*tl_Tstar(i))/TairK(i)-   &
#ifdef TL_IOMS
     &                 upvel
#endif
            Hlw=rhoAir(i)*Hlv(i,j)*upvel*Q(i)
            tl_Hlw=rhoAir(i)*Q(i)*(tl_Hlv(i,j)*upvel+                   &
     &                         Hlv(i,j)*tl_upvel)-                      &
#ifdef TL_IOMS
     &                         Hlw
#endif

            LHeat(i,j)=(Hl+Hlw)
            tl_LHeat(i,j)=(tl_Hl+tl_Hlw)
# ifdef MASKING
            LHeat(i,j)=LHeat(i,j)*rmask(i,j)
            tl_LHeat(i,j)=tl_LHeat(i,j)*rmask(i,j)
# endif
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
            Taur=0.85_r8*rain(i,j)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
            cff=rhoAir(i)*Cd*Wspeed
            tl_cff=rhoAir(i)*(tl_Cd*Wspeed+Cd*tl_Wspeed)-               &
#ifdef TL_IOMS
     &                         cff
#endif
            Taux(i,j)=(cff*Uwind(i,j)+Taur*SIGN(1.0_r8,Uwind(i,j)))
            tl_Taux(i,j)=tl_cff*Uwind(i,j)+                             &
#ifdef TL_IOMS
     &                      Taur*SIGN(1.0_r8,Uwind(i,j))
#endif
# ifdef MASKING
            Taux(i,j)=Taux(i,j)*rmask(i,j)
            tl_Taux(i,j)=tl_Taux(i,j)*rmask(i,j)
# endif
            Tauy(i,j)=(cff*Vwind(i,j)+Taur*SIGN(1.0_r8,Vwind(i,j)))
            tl_Tauy(i,j)=tl_cff*Vwind(i,j)+                             &
#ifdef TL_IOMS
     &                      Taur*SIGN(1.0_r8,Vwind(i,j))
#endif
# ifdef MASKING
            Tauy(i,j)=Tauy(i,j)*rmask(i,j)
            tl_Tauy(i,j)=tl_Tauy(i,j)*rmask(i,j)
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
!  to (psu m/s) for stflx(isalt) in "set_vbc.F". The E-P value is
!  saved in variable EminusP for I/O purposes.
!
      Hscale=1.0_r8/(rho0*Cp)
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          lrflx(i,j)=LRad(i,j)*Hscale
          tl_lrflx(i,j)=tl_LRad(i,j)*Hscale
          lhflx(i,j)=-LHeat(i,j)*Hscale
          tl_lhflx(i,j)=-tl_LHeat(i,j)*Hscale
          shflx(i,j)=-SHeat(i,j)*Hscale
          tl_shflx(i,j)=-tl_SHeat(i,j)*Hscale
          stflx(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                      &
     &                      lhflx(i,j)+shflx(i,j))
          tl_stflx(i,j,itemp)=(tl_lrflx(i,j)+                           &
     &                      tl_lhflx(i,j)+tl_shflx(i,j))+               &
# ifdef TL_IOMS
     &                       srflx(i,j)
# endif
# ifdef EMINUSP
          evap(i,j)=LHeat(i,j)/Hlv(i,j)
          tl_evap(i,j)=tl_LHeat(i,j)/Hlv(i,j)-                          &
     &                     tl_Hlv(i,j)*evap(i,j)/Hlv(i,j)+              &
#  ifdef TL_IOMS
     &                     evap(i,j)
#  endif
          stflx(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
          tl_stflx(i,j,isalt)=cff*tl_evap(i,j)-                         &
#  ifdef TL_IOMS
     &                       cff*rain(i,j)
#  endif
# endif
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          tl_stflx(i,j,itemp)=tl_stflx(i,j,itemp)*rmask(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask(i,j)
          tl_evap(i,j)=tl_evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
          tl_stflx(i,j,isalt)=tl_stflx(i,j,isalt)*rmask(i,j)
#   endif
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
          tl_sustr(i,j)=cff*(tl_Taux(i-1,j)+tl_Taux(i,j))
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
          tl_sustr(i,j)=tl_sustr(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
          tl_svstr(i,j)=cff*(tl_Tauy(i,j-1)+tl_Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
          tl_svstr(i,j)=tl_svstr(i,j)*vmask(i,j)
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
     &                          tl_lrflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_stflx(:,:,itemp))
# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_stflx(:,:,isalt))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_evap)
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_svstr)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iRPM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iRPM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_lrflx, tl_lhflx, tl_shflx,                 &
     &                    tl_stflx(:,:,itemp))
#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap,                                         &
     &                    stflx(:,:,isalt))
      CALL mp_exchange2d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_evap,                                      &
     &                    tl_stflx(:,:,isalt))
#  endif
      CALL mp_exchange2d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
      CALL mp_exchange2d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_sustr, tl_svstr)
# endif

      RETURN
      END SUBROUTINE rp_bulk_flux_tile
#endif
      END MODULE rp_bulk_flux_mod
