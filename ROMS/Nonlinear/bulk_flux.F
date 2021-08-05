#include "cppdefs.h"
      MODULE bulk_flux_mod
#ifdef BULK_FLUXES
!
!svn $Id: bulk_flux.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!  Modified by Jim Edson to correct specific humidities.               !
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
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
# ifdef EMINUSP
     &                     FORCES(ng) % evap,                           &
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
     &                           rain, lhflx, lrflx, shflx, srflx,      &
# ifdef EMINUSP
     &                           evap,                                  &
# endif
     &                           stflux, sustr, svstr)
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
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  endif
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

      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflux(LBi:,LBj:,:)

#  ifdef EMINUSP
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
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
#  endif
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

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflux(LBi:UBi,LBj:UBj,NT(ng))

#  ifdef EMINUSP
      real(r8), intent(out) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
      integer, parameter :: IterMax = 3

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Cd, Hl, Hlw, Hscale, Hscale2, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL

      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
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

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
#ifdef WIND_MINUS_CURRENT
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
#else
!
!  Load wind components to local arrays.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Uair(i,j)=Uwind(i,j)
          Vair(i,j)=Vwind(i,j)
        END DO
      END DO
#endif
!
!  Compute Atmosphere-ocean fluxes using a bulk flux parameterization.
!
      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
      twopi_inv=0.5_r8/pi
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uair(i,j)*Uair(i,j)+Vair(i,j)*Vair(i,j))+eps
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
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff+eps)) !Spec hum (kg/kg)
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
#ifdef REGCM_COUPLING
          Hlv(i,j)=2.5104E+6_r8 ! to be consistent with RegCM
#else
          Hlv(i,j)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
#endif
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
# if defined DRAGLIM_DAVIS
            ZoW(i)=MAX(ZoW(i),1.27E-7_r8) ! Davis et al. (2008)
            ZoW(i)=MIN(ZoW(i),2.85E-3_r8)
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
          diffh=0.02411_r8*(1.0_r8+TairC(i)*                            &
     &                      (3.309E-3_r8-1.44E-6_r8*TairC(i)))/         &
     &          (rhoAir(i)*blk_Cpa+eps)
          cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
          wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/        &
     &                                     (blk_Cpa*diffh))
          Hsr=rain(i,j)*wet_bulb*blk_Cpw*                               &
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
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
          Taur=0.85_r8*rain(i,j)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
          cff=rhoAir(i)*Cd*Wspeed
          Taux(i,j)=(cff*Uair(i,j)+Taur*SIGN(1.0_r8,Uair(i,j)))
# ifdef MASKING
          Taux(i,j)=Taux(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          Taux(i,j)=Taux(i,j)*rmask_wet(i,j)
# endif
          Tauy(i,j)=(cff*Vair(i,j)+Taur*SIGN(1.0_r8,Vair(i,j)))
# ifdef MASKING
          Tauy(i,j)=Tauy(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          Tauy(i,j)=Tauy(i,j)*rmask_wet(i,j)
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
!  Note that stflux(:,:,isalt) is the E-P flux. The fresh water flux
!  is positive out of the ocean and the salt flux is positive into the
!  ocean. It is  multiplied by surface salinity when computing state
!  variable stflx(:,:,isalt) in "set_vbc.F".
!
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          lrflx(i,j)=LRad(i,j)*Hscale2
          lhflx(i,j)=-LHeat(i,j)*Hscale2
          shflx(i,j)=-SHeat(i,j)*Hscale2
          stflux(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                     &
     &                       lhflx(i,j)+shflx(i,j))
# ifdef MASKING
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask_wet(i,j)
# endif
# ifdef EMINUSP
          evap(i,j)=LHeat(i,j)/(Hlv(i,j)+eps)
#  ifdef MASKING
          evap(i,j)=evap(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          evap(i,j)=evap(i,j)*rmask_wet(i,j)
#  endif
          stflux(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
#  ifdef MASKING
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask_wet(i,j)
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
# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflux(:,:,isalt))
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflux(:,:,itemp))
#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap,                                         &
     &                    stflux(:,:,isalt))
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
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
#endif
      END MODULE bulk_flux_mod
