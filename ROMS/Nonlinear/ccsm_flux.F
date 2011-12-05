#include "cppdefs.h"
      MODULE ccsm_flux_mod
#if defined CCSM_FLUXES && defined BULK_FLUXES
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!       CCSM's flux_mod.F90 by B. Kauffman                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC  :: ccsm_flux

contains
!
!***********************************************************************
      SUBROUTINE ccsm_flux (ng, tile)
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
      CALL ccsm_flux_tile (ng, tile,                                    &
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
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % beta,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
# if defined CLOUDS
     &                     FORCES(ng) % cloud,                          &
# endif
# ifdef COARE_TAYLOR_YELLAND
     &                     FORCES(ng) % Hwave,                          &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
     &                     FORCES(ng) % Pwave_top,                      &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
     &                     FORCES(ng) % Lwave,                          &
# endif
# ifdef RUNOFF
     &                     FORCES(ng) % runoff,                         &
# endif
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
     &                     FORCES(ng) % stflx,                          &
# if defined SHORTWAVE && defined ALBEDO
     &                     FORCES(ng) % cawdir,                         &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                     GRID(ng) % latr,                             &
# endif
# ifdef ICE_MODEL
     &                     ICE(ng) % ai,                                &
     &                     ICE(ng) % hi,                                &
#  ifdef ICE_THERMO
     &                     FORCES(ng) % qai_n,                          &
     &                     ICE(ng) % hsn,                               &
     &                     ICE(ng) % tis,                               &
     &                     ICE(ng) % sfwat,                             &
     &                     ICE(ng) % coef_ice_heat,                     &
     &                     ICE(ng) % rhs_ice_heat,                      &
     &                     FORCES(ng) % snow_n,                         &
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
      END SUBROUTINE ccsm_flux
!
!***********************************************************************
      SUBROUTINE ccsm_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liold, linew,                          &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
     &                           alpha, beta, rho, t, u, v,             &
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
# if defined CLOUDS
     &                           cloud,                                 &
# endif
# ifdef COARE_TAYLOR_YELLAND
     &                           Hwave,                                 &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
     &                           Pwave_top,                             &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
     &                           Lwave,                                 &
# endif
# ifdef RUNOFF
     &                           runoff,                                &
# endif
     &                           rain, lhflx, lrflx, shflx,             &
     &                           srflx, stflx,                          &
# if defined SHORTWAVE && defined ALBEDO
     &                           cawdir,                                &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                           latr,                                  &
# endif
# ifdef ICE_MODEL
     &                           ai, hi,                                &
#  ifdef ICE_THERMO
     &                           qai_n,                                 &
     &                           hsn, tis, sfwat,                       &
     &                           coef_ice_heat, rhs_ice_heat, snow_n,   &
#  endif
     &                           sustr_aw, svstr_aw, qao_n,             &
     &                           tau_aix_n, tau_aiy_n,                  &
# endif
# ifdef EMINUSP
     &                           EminusP, evap,                         &
# endif
     &                           sustr, svstr)
!***********************************************************************
!     call srfflx_ao( nloc    ,   z   ,   u   ,   v   ,ptem   , &
!     &              shum     ,dens   ,uocn   ,vocn   , &
!     &              tocn     ,mask   , shflx   , lhflx   ,lwup   , &
!     &              evap     ,taux   , tauy  ,tref   ,qref   , &
!     &              duu10n                                     )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
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
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
#  ifdef COARE_TAYLOR_YELLAND
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
#  endif
#  if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: rain(LBi:,LBj:)
#  ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
#  if defined SHORTWAVE && defined ALBEDO
      real(r8), intent(inout) :: cawdir(LBi:,LBj:)
#  endif
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(inout) :: latr(LBi:,LBj:)
#  endif

#  ifdef ICE_MODEL
      real(r8), intent(out) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(out) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(out) :: qao_n(LBi:,LBj:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qai_n(LBi:,LBj:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
      real(r8), intent(in) :: sfwat(LBi:,LBj:,:)
      real(r8), intent(out) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: rhs_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: snow_n(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_aiy_n(LBi:,LBj:)
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
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  if defined CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
# ifdef COARE_TAYLOR_YELLAND
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
#  ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# if defined SHORTWAVE && defined ALBEDO
      real(r8), intent(inout) :: cawdir(LBi:UBi,LBj:UBj)
#  endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(inout) :: latr(LBi:UBi,LBj:UBj)
#  endif

#  ifdef ICE_MODEL
      real(r8), intent(out) :: sustr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qao_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qai_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sfwat(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: coef_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhs_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: snow_n(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: tau_aix_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tau_aiy_n(LBi:UBi,LBj:UBj)
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
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: Iter, i, j, k
      integer :: IterMax

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Ce, Ch, Hl, Hlw, Hscale, Hscale2, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur, cff1, cff2
      real(r8) :: Wspeed, ZQoL, ZToL, cff
# ifdef ICE_BULK_FLUXES
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
      real(r8), parameter :: alb_w=0.1_r8
      real(r8) :: albs, qswi, qlwi, qlh_i, qsh_i
      real(r8) :: le_i, dq_i,fqlat1, sfc_temp, slp, Qsati
      real(r8) :: vap_p_i
#  ifdef ICE_ALB_EC92
      real(r8) :: albi, albsn, thki_n, thksn_n
#  endif
# endif
# if defined SHORTWAVE && defined ALBEDO
      real(r8) :: cawdif
      real(r8), parameter :: albw_d=0.065_r8
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8) :: caw_curve
# endif
# if defined ICE_MODEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
# endif

! ------------------------------------------------------------------
! !DESCRIPTION:
!     wrapper to atm/ocn flux calculation
!
! !REMARKS:
!     All data must be on the ocean domain (note: a domain includes a
!     particular decomposition).
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - first version
!
! ------------------------------------------------------------------
!
! !IROUTINE: srfflx_ao -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - brought in from cpl5.
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!
! !INTERFACE: ------------------------------------------------------------------

!SUBROUTINE srfflx_ao_tile(blk_ZX  ,Uwind  ,Vwind  ,Tair ,   &
!          &         Hair  ,rhoAir  ,u    ,v    ,   &
!          &         t    ,rmask  ,shflx   ,lhflx   ,lwup  ,   &
!          &         evap  ,taux  ,tauy  ,tref  ,qref  ,   &
!          &         duu10n                                )

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
!  integer(IN),intent(in) :: rmask(imax) ! ocn domain mask 0
!  real(R8)   ,intent(in) :: blk_ZX (imax) ! atm level height      (m)
!  real(R8)   ,intent(in) :: Uwind (imax) ! atm u wind            (m/s)
!  real(R8)   ,intent(in) :: Vwind (imax) ! atm v wind            (m/s)
!  real(R8)   ,intent(in) :: Tair(imax) ! atm potential T       (K)
!  real(R8)   ,intent(in) :: Hair (imax) ! atm specific humidity (kg/kg)
!  real(R8)   ,intent(in) :: rhoAir (imax) ! atm air density       (kg/m^3)
!  real(R8)   ,intent(in) :: u   (imax) ! ocn u-velocity        (m/s)
!  real(R8)   ,intent(in) :: v   (imax) ! ocn v-velocity        (m/s)
!  real(R8)   ,intent(in) :: t   (imax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
!  real(R8),intent(out)  ::  shflx  (imax) ! heat flux: sensible    (W/m^2)
!  real(R8),intent(out)  ::  lhflx  (imax) ! heat flux: latent      (W/m^2)
!  real(R8),intent(out)  ::  lwup (imax) ! heat flux: lw upward   (W/m^2)
!  real(R8),intent(out)  ::  evap (imax) ! water flux: evap  ((kg/s)/m^2)
!  real(R8),intent(out)  ::  taux (imax) ! surface stress, zonal      (N)
!  real(R8),intent(out)  ::  tauy (imax) ! surface stress, maridional (N)
!  real(R8),intent(out)  ::  tref (imax) ! diag:  2m ref height T     (K)
!  real(R8),intent(out)  ::  qref (imax) ! diag:  2m ref humidity (kg/kg)
!  real(R8),intent(out)  :: duu10n(imax) ! diag: 10m wind speed squared (m/s)^2

!EOP

   !--- local constants --------------------------------
      real(R8),parameter :: umin  =  0.5    ! minimum wind speed       (m/s)
      real(R8),parameter :: zref  = 10.0    ! reference height           (m)
      real(R8),parameter :: ztref =  2.0    ! reference height for air T (m)

   !--- local variables --------------------------------
      real(r8)    :: lwup   ! upward longwave radiation
      real(r8)    :: vmag   ! surface wind magnitude   (m/s)
      real(r8)    :: thvbot ! virtual temperature      (K)
      real(r8)    :: ssq    ! sea surface humidity     (kg/kg)
      real(r8)    :: delt   ! potential T difference   (K)
      real(r8)    :: delq   ! humidity difference      (kg/kg)
      real(r8)    :: stable ! stability factor
      real(r8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
      real(r8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
      real(r8)    :: ren    ! sqrt of neutral exchange coeff (water)
      real(r8)    :: rd     ! sqrt of exchange coefficient (momentum)
      real(r8)    :: re     ! sqrt of exchange coefficient (water)
      real(r8)    :: ustar  ! ustar
      real(r8)    :: qstar  ! qstar
      real(r8)    :: tstar  ! tstar
      real(r8)    :: hol    ! H (at blk_ZX) over L
      real(r8)    :: xsq    ! ?
      real(r8)    :: xqq    ! ?
      real(r8)    :: psimh  ! stability function at blk_ZX (momentum)
      real(r8)    :: psixh  ! stability function at blk_ZX (heat and water)
      real(r8)    :: psix2  ! stability function at ztref reference height
      real(r8)    :: alz    ! ln(blk_ZX/zref)
      real(r8)    :: al2    ! ln(zref/ztref)
      real(r8)    :: u10n   ! 10m neutral wind
      real(r8)    :: tau    ! stress at blk_ZX
      real(r8)    :: cpair     ! specific heat of moist air
      real(r8)    :: bn     ! exchange coef funct for interpolation
      real(r8)    :: bh     ! exchange coef funct for interpolation
      real(r8)    :: fac    ! vertical interpolation factor
      real(r8)    :: Hlv    ! Latent heat of vaporization at sea surface
      real(r8), dimension(IminS:ImaxS) :: rhoAir
      real(r8), dimension(IminS:ImaxS) :: TairK
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
#endif

   !--- local functions --------------------------------
      real(r8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
      real(r8)    :: cdn    ! function: neutral drag coeff at 10m
      real(r8)    :: psimhu ! function: unstable part of psimh
      real(r8)    :: psixhu ! function: unstable part of psimx
      real(r8)    :: tsk    ! water temperature (K)
      real(r8)    :: Tk     ! dummy arg ~ air temperature (K)
      real(r8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
      real(r8)    :: xd     ! dummy arg ~ ?

!  Stefan-Boltzmann constant ~ W/m^2/K^4
      real(r8),parameter :: shr_const_stebol  = 5.67e-8_r8
! Boltzmann's constant ~ J/K/molecule
      real(r8),parameter :: shr_const_boltz   = 1.38065e-23_r8
! Avogadro's number ~ molecules/kmole
      real(r8),parameter :: shr_const_avogad  = 6.02214e26_r8
! Universal gas constant ~ J/K/kmole
      real(r8),parameter :: shr_const_rgas =                            &
     &         shr_const_avogad*shr_const_boltz
! molecular weight dry air ~ kg/kmole
      real(r8),parameter :: shr_const_mwdair  = 28.966_r8
! molecular weight water vapor
      real(r8),parameter :: shr_const_mwwv = 18.016_r8
! Water vapor gas constant ~ J/K/kg
      real(r8),parameter :: shr_const_rwv =                             &
     &         shr_const_rgas/shr_const_mwwv
! Dry air gas constant     ~ J/K/kg
      real(r8),parameter :: shr_const_rdair   =                         &
     &         shr_const_rgas/shr_const_mwdair
      real(r8),parameter :: shr_const_zvir    =                         &
     &         (shr_const_rwv/shr_const_rdair)-1.0_r8
! Von Karman constant
      real(r8),parameter :: shr_const_karman  = 0.4_r8
! specific heat of dry air   ~ J/kg/K
      real(r8),parameter :: shr_const_cpdair  = 1.00464e3_r8
! specific heat of water vap ~ J/kg/K
      real(r8),parameter :: shr_const_cpwv    = 1.810e3_R8
! CPWV/CPDAIR - 1.0
      real(r8),parameter :: shr_const_cpvir   =                          &
     &         (shr_const_cpwv/shr_const_cpdair)-1.0_r8
! latent heat of evaporation ~ J/kg
      real(r8),parameter :: shr_const_latvap  = 2.501e6_r8


      qsat(Tk)   = 640380.0 / exp(5107.4/Tk)
      cdn(Umps)  = 0.0027 / Umps + 0.000142 + 0.0000764 * Umps
      psimhu(xd) = log((1.0+xd*(2.0+xd))*(1.0+xd*xd)/8.0) - 2.0*atan(xd) + 1.571
      psixhu(xd) = 2.0 * log((1.0 + xd*xd)/2.0)

#include "set_bounds.h"
!-------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES:
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + shflx + lhflx
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------

      al2 = log(zref/ztref)

      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
! Note that this loop needs to be cleaned of all global arrays for
! OpenMP.
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          !--- compute some needed quantities ---
          vmag   = max(umin, sqrt( (Uwind(i,j)-u(i,j,N(ng),nrhs))**2 +  &
     &                             (Vwind(i,j)-v(i,j,N(ng),nrhs))**2) )
          TairK(i) = Tair(i,j) + 273.16_r8
          thvbot = TairK(i) * (1.0 + shr_const_zvir * Hair(i,j)) ! virtual temp (K)
          rhoAir(i) = Pair(i,j)*100.0_r8/(blk_Rgas*TairK(i)*                  &
     &                    (1.0_r8+0.61_r8*Hair(i,j))) !  Moist air density (kg/m3).
          ssq    = 0.98_r8 * qsat(TairK(i)) / rhoAir(i)    ! sea surf hum (kg/kg)
          delt   = Tair(i,j) - t(i,j,N(ng),nrhs,itemp)                  ! pot temp diff (K)
          delq   = Hair(i,j) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(blk_ZW(ng)/zref)
          cpair     = shr_const_cpdair*(1.0 + shr_const_cpvir*ssq)

          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------

          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_r8 + sign(0.5_r8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8
          ren    = 0.0346_r8

          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq

          !--- compute stability & evaluate all stability functions ---
          hol  = shr_const_karman*g*blk_ZT(ng)*                         &
     &          (tstar/thvbot+qstar/(1.0_r8/shr_const_zvir+Hair(i,j)))/ &
     &                     ustar**2
          hol  = sign( min(abs(hol),10.0_r8), hol )
          stable = 0.5_r8 + sign(0.5_r8 , hol)
          xsq    = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)

          !--- shift wind speed using old coefficient ---
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd / rdn

          !--- update transfer coeffs at 10m and neutral stability ---
          rdn = sqrt(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8-stable)*0.0327_r8 + stable * 0.018_r8

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))

          !--- update ustar, tstar, qstar using updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !------------------------------------------------------------
          ! iterate to converge on Z/L, ustar, tstar and qstar
          !------------------------------------------------------------

          !--- compute stability & evaluate all stability functions ---
          hol  = shr_const_karman*g*blk_ZQ(ng)*                         &
     &            (tstar/thvbot+qstar/(1.0/shr_const_zvir+Hair(i,j)))/  &
     &                         ustar**2
          hol  = sign( min(abs(hol),10.0_r8), hol )
          stable = 0.5_r8 + sign(0.5_r8 , hol)
          xsq    = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)

          !--- shift wind speed using old coeffs ---
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd/rdn

          !--- update transfer coeffs at 10m and neutral stability ---
          rdn = sqrt(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8 - stable)*0.0327_r8 + stable * 0.018_r8

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))

          !--- update ustar, tstar, qstar using updated, shifted coeffs ---
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------

          tau = rhoAir(i) * ustar * ustar
          SRad(i,j) = srflx(i,j) * Hscale

          !--- momentum flux ---
          Taux(i,j) = tau * (Uwind(i,j)-u(i,j,N(ng),nrhs)) / vmag
          Tauy(i,j) = tau * (Vwind(i,j)-v(i,j,N(ng),nrhs)) / vmag
# ifdef ICE_MODEL
          Taux_Ice(i,j) = rhoAir(i)*cdai(ng)*Uwind(i,j)*vmag
          Tauy_Ice(i,j) = rhoAir(i)*cdai(ng)*Vwind(i,j)*vmag

#  ifdef ICE_THERMO
!  Correct ocean net short-wave radiation for ice cover and convert to
!  kinematic units

          Ch=ustar*tstar/(vmag*delt+0.0098_r8*blk_ZT(ng))
          Ce=ustar*qstar/(vmag*delq+eps)

          Coef_Ice(i,j) = 0._r8
          RHS_Ice(i,j) = 0._r8

! Sensible heat flux over ice

          IF (ai(i,j,liold).gt.min_a(ng)) THEN
            sfc_temp = tis(i,j) + 273.16_r8
          ELSE
            sfc_temp = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          ENDIF
          qsh_i = rhoAir(i)*blk_Cpa*Ch*vmag*                            &
     &            (TairK(i)+0.0098*blk_ZT(ng)-sfc_temp)

          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &          rhoAir(i)*blk_Cpa*Ch*vmag
          RHS_Ice(i,j) = RHS_Ice(i,j) +                                 &
     &        rhoAir(i)*blk_Cpa*Ch*vmag*(TairK(i)+0.0098*blk_ZT(ng))

! Latent heat flux over ice
          le_i=2.834E+6_r8

! Qsati is saturation specific humidity over the ice from
! Parkinson and Washington (1979)
          slp = Pair(i,j)*100._r8 !Convert back to Pascal from millibars
          cff = 611._r8*                                                &
     &      10._r8**(9.5_r8*(sfc_temp-273.16_r8)/(sfc_temp-7.66_r8))
          Qsati=eps*cff/(slp-(1._r8-.62197_r8)*cff)
          dq_i=Hair(i,j)-Qsati
          fqlat1=rhoAir(i)*Ce*vmag
          qlh_i = fqlat1*le_i*dq_i
          RHS_Ice(i,j) = RHS_Ice(i,j) + qlh_i
!
! Calculate the ice/snow albedo
#   ifdef ICE_ALB_EC92
! Ice and snow albedo is calculated from Ebert and Curry,1992b
          albi=0._r8
          albsn=0._r8
          IF (ai(i,j,liold) .gt. min_a(ng)) THEN
            thki_n = hi(i,j,liold)/(ai(i,j,liold)+0.01_r8)
            thksn_n = hsn(i,j,liold)/(ai(i,j,liold)+0.01_r8)
            albi=0.082409_r8*LOG(thki_n)+0.485472_r8
            IF (thki_n.GE.1._r8) albi=0.07616_r8*thki_n+0.414492_r8
            IF (thki_n.GE.2._r8) albi=0.561632_r8
!approximated values for albsn(depends on COSZ, but small variation)
            albsn=0.83_r8
            IF (sfwat_n.GT.0._r8) THEN
              IF (thksn_n.GE.0.1_r8) THEN
                albsn=0.701009_r8
              ELSE
                albsn=albi+(0.701009_r8-albi)*thksn_n/0.1
              ENDIF
            ENDIF
            albs=albi
            IF (hsn(i,j,liold).GT.0._r8) albs=albsn
!            IF (sfwat_n.GT.0._r8) albs=0.100737_r8                      &
!     &         +0.518_r8*EXP(-8.1_r8 *sfwat_n-0.47_r8)                  &
!     &         +0.341_r8*EXP(-31.8_r8*sfwat_n-0.94_r8)                  &
!     &         +0.131_r8*EXP(-2.6_r8 *sfwat_n-3.82_r8)

            ENDIF
          ELSE
            albs=alb_w
          ENDIF
#   else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (ai(i,j,liold) .gt. min_a(ng)) THEN
            IF (hsn(i,j,liold).gt.0._r8) THEN
              IF (sfc_temp-273.16_r8.gt. -1._r8) THEN
                albs = cff1*(sfc_temp-272.16_r8)+alb_s_dry
              ELSE
                albs = alb_s_dry
              ENDIF
            ELSE
              IF (sfc_temp-273.16_r8.gt. -1._r8) THEN
                albs = cff2*(sfc_temp-272.16_r8)+alb_i_dry
              ELSE
                albs = alb_i_dry
              ENDIF
            ENDIF
          ELSE
            albs=alb_w
          ENDIF
#   endif
!
! Short-wave radiation to ice
          qswi = (1._r8-albs)*SRad(i,j)
! Already in W/m**2
          RHS_Ice(i,j) = RHS_Ice(i,j) + qswi

! Net long-wave radiation from ice
          cff = Hair(i,j)/(1.0_r8+Hair(i,j))
          vap_p_i = slp*cff/(0.62197 + cff)
# ifdef LONGWAVE
          qlwi = emmiss*StefBo*(TairK(i)**4)*                           &
     &         (0.39_r8-0.005_r8*SQRT(vap_p_i))*                        &
     &                (1.0_r8-0.65_r8*cloud(i,j)) +                     &
     &          4._r8*StefBo*emmiss*TairK(i)**3*(sfc_temp-TairK(i))

          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &           4._r8*StefBo*emmiss*TairK(i)**3
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &         emmiss*StefBo*(TairK(i)**4)*(                            &
     &    (0.39_r8-0.005_r8*SQRT(vap_p_i))*(1.0_r8-0.65_r8*cloud(i,j))  &
     &    - 4.0_r8)
# else
          qlwi = lrflx(i,j)*Hscale
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &           lrflx(i,j)*Hscale -                                    &
     &     4._r8*StefBo*emmiss*TairK(i)**3*(sfc_temp-2._r8*TairK(i))

          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &            4._r8*StefBo*emmiss*TairK(i)**3
# endif
! Correct for Kelvin to Celsius in RHS
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &            Coef_Ice(i,j)*273.16_r8

! Net heat flux from ice to atmosphere
          Qai(i,j) =-qsh_i - qlh_i - qswi + qlwi
#  endif
# endif
          !--- heat flux ---
          SHeat (i,j) = cpair * tau * tstar / ustar
          LHeat (i,j) = shr_const_latvap * tau * qstar / ustar
        END DO
      END DO

      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          shflx (i,j) = SHeat(i,j)
          lhflx (i,j) = LHeat(i,j)
          tsk = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          lwup = -shr_const_stebol * tsk**4

# if defined SHORTWAVE && defined ALBEDO
          cawdif = 1._r8-albw_d
          srflx(i,j) = srflx(i,j)*(cawdir(i,j)*                         &
     &                 (1._r8-cloud(i,j))+cawdif*cloud(i,j))
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
          caw_curve = 1._r8-(0.069_r8 - 0.011_r8*                       &
     &                        cos(2*deg2rad*latr(i,j)))
          srflx(i,j) = srflx(i,j)*caw_curve
# endif
          lrflx(i,j) = lwup*Hscale2 + lrflx(i,j)
          lhflx(i,j) = lhflx(i,j)*Hscale2
          shflx(i,j) = shflx(i,j)*Hscale2
          stflx(i,j,itemp) = (srflx(i,j)+lrflx(i,j)+                    &
     &                      lhflx(i,j)+shflx(i,j))
# ifdef EMINUSP
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
          Hlv = (2.501_r8-0.00237_r8*t(i,j,N(ng),nrhs,itemp))*1.0e+6_r8
          evap(i,j) = -LHeat(i,j)/Hlv
          stflx(i,j,isalt) = cff*(evap(i,j)-rain(i,j))
#  ifdef ICE_MODEL
          evap(i,j) = (1.0_r8-ai(i,j,linew))*evap(i,j)
          stflx(i,j,isalt) = (1.0_r8-ai(i,j,linew))*stflx(i,j,isalt)
#  endif
# endif
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#   endif
# endif
# ifdef RUNOFF
          stflx(i,j,isalt) = stflx(i,j,isalt) -                         &
     &                        cff*runoff(i,j)
#  ifdef MASKING
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#  endif
# endif

# if defined ICE_MODEL
! Limit shortwave radiation to be used in computing penetrative
! radiation by presence of ice and biology

          ice_thick(i,j) = hi(i,j,linew)/(ai(i,j,linew)+eps)
          snow_thick(i,j) = hsn(i,j,linew)/(ai(i,j,linew)+eps)

          IF (ice_thick(i,j).le.0.1_r8) THEN
            cff1=ice_thick(i,j)*5._r8
          ELSE
            cff1=0.1_r8*5._r8+(ice_thick(i,j)-0.1_r8)*1._r8
          ENDIF

          cff1= cff1+snow_thick(i,j)*20._r8
          srflx(i,j)=(1-ai(i,j,linew))*srflx(i,j)                       &
     &            +ai(i,j,linew)*srflx(i,j)*exp(-cff1)

          qao_n(i,j) = -stflx(i,j,itemp)*Hscale

! Net precipitation - evaporation in m/s
          IF (Tair(i,j).ge.0.0) THEN
            snow_n(i,j) = 0.0_r8
          ELSE
            snow_n(i,j) = rain(i,j)/rhosnow_dry(ng)
          ENDIF
          coef_ice_heat(i,j) = Coef_Ice(i,j)
          rhs_ice_heat(i,j) = RHS_Ice(i,j)
          qai_n(i,j) = Qai(i,j)

# endif

!--- water flux ---
!         evap(i,j) = lhflx(i,j)/shr_const_latvap

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
!         hol = hol*ztref/blk_ZQ(ng)
!         xsq = max( 1.0, sqrt(abs(1.0-16.0*hol)) )
!         xqq = sqrt(xsq)
!         psix2   = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
!         fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
!         tref(i,j) = TairK - delt*fac
!         tref(i,j) = tref(i,j) - 0.01*ztref   ! pot temp to temp correction
!         fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
!         qref(i,j) =  Hair(i,j) - delq*fac

!         duu10n(i,j) = u10n*u10n ! 10m wind speed squared

# ifdef EMINUSP
          EminusP(i,j)=stflx(i,j,isalt)
# endif
        END DO
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      Hscale = 1./rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=Taux(i,j)*Hscale
# ifdef MASKING
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
          svstr(i,j)=Tauy(i,j)*Hscale
# ifdef MASKING
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
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        lrflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        lhflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        shflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,itemp))
#  ifdef ICE_MODEL
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        srflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        qao_n)
#   ifdef ICE_THERMO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    qai_n)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    coef_ice_heat)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    rhs_ice_heat)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    snow_n)
#   endif
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tau_aix_n)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tau_aiy_n)
#  endif
#  ifdef EMINUSP
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        evap)
       CALL exchange_r2d_tile (ng, tile,                                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        EminusP)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,isalt))
#  endif
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr)
#  ifdef ICE_MODEL
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr_aw)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr_aw)
#  endif
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflx(:,:,itemp))
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    srflx, qao_n, Tau_aix_n, Tau_aiy_n )
#   ifdef ICE_THERMO
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    qai_n, coef_ice_heat, rhs_ice_heat, snow_n )
#   endif
#  endif
#  if defined EMINUSP || (defined ICE_MODEL && defined ICE_THERMO)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    evap, EminusP,                                &
     &                    stflx(:,:,isalt))
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sustr, svstr)
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sustr_aw, svstr_aw)
#  endif
# endif

      RETURN
      END SUBROUTINE ccsm_flux_tile
#endif

      END module ccsm_flux_mod
