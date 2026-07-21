#include "cppdefs.h"
      MODULE marsh_vert_growth_mod

#if defined MARSH_VERT_GROWTH
!
!git $Id$
!=============================================== Tarandeep S. Kalra ===!
!  Copyright (c) 2002-2026 The ROMS Group                 Joel Carr    !
!    Licensed under a MIT/X style license             Neil K. Ganju    !
!    See License_ROMS.md                       Alfredo Aretxabaleta    !
!==================================================  John C. Warner ====
!                                                                      !
!  It Calculates vertical growth of marsh through biomass production.  !
!  Biomass production can lead to mudflats being converted to marsh    !
!  cells and evolve marsh stem/canopy properties                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    D'Alpaos, A.,  Lanzoni, S.,  Marani, M., Rinaldo, A., 2006:       !
!      Landscape evolution in tidal embayments: Modeling the interplay !
!      of erosion, sedimentation, and vegetation dynamics, J. Geophys. !
!      Res., 112, doi: 10.1029/2006JF000537.                           !
!      [Used for marsh vegetation parameters]                          !
!                                                                      !
!    Morris, J. T., P. V. Sundareshwar, C. T. Nietch, B. Kjerfve, and  !
!      D. R. Cahoon, 2002: Responses of coastal wetlands to rising     !
!      sea level, Ecology, 83 (10), 2869-2877.                         !
!      [Used for marsh vertical growth and parabolic profile]          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_sedbed
      USE mod_sediment
      USE mod_vegetation
!
      USE bc_2d_mod,       ONLY : bc_r2d_tile
      USE bc_3d_mod,       ONLY : bc_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d,         &
     &                            mp_exchange4d
# endif
!
      implicit none
!
      PUBLIC  :: marsh_vert_growth
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE marsh_vert_growth (ng,tile)
!***********************************************************************
!
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 44, __LINE__, MyFile)
# endif
      CALL marsh_vert_growth_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nnew(ng),                            &
     &                             GRID(ng) % h,                        &
     &                             OCEAN(ng)  % zeta,                   &
     &                             SEDBED(ng) % bed,                    &
     &                             SEDBED(ng) % bed_mass,               &
     &                             VEG(ng) % marsh_mask,                &
!!# if defined MARSH_TIDAL_RANGE
     &                             VEG(ng) % marsh_high_water,          &
!!#  if defined MARSH_MCKEE_FORMULATION
     &                             VEG(ng) % marsh_low_water,           &
!!#  endif
!!   &                             VEG(ng) % marsh_tidal_range,         &
!!# endif
# if defined MARSH_BIOMASS_VEG
     &                             VEG(ng) % plant,                     &
# endif
     &                             VEG(ng) % marsh_biomass_peak,        &
     &                             VEG(ng) % marsh_vert_rate,           &
     &                             VEG(ng) % marsh_accret)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 44, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE marsh_vert_growth

!***********************************************************************
      SUBROUTINE marsh_vert_growth_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nnew,                          &
     &                                   h,                             &
     &                                   zeta,                          &
     &                                   bed,                           &
     &                                   bed_mass,                      &
     &                                   marsh_mask,                    &
!!# if defined MARSH_TIDAL_RANGE
     &                                   marsh_high_water,              &
!!#  if defined MARSH_MCKEE_FORMULATION
     &                                   marsh_low_water,               &
!!#  endif
!!   &                                   marsh_tidal_range,             &
!!# endif
# if defined MARSH_BIOMASS_VEG
     &                                   plant,                         &
# endif
     &                                   marsh_biomass_peak,            &
     &                                   marsh_vert_rate,               &
     &                                   marsh_accret)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: marsh_mask(LBi:,LBj:)
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(in)    :: zeta(LBi:,LBj:,:)
      real(r8), intent(in)    :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(in)    :: marsh_high_water(LBi:,LBj:)
!!#  if defined MARSH_MCKEE_FORMULATION
      real(r8), intent(in)    :: marsh_low_water(LBi:,LBj:)
!!#  endif
!!    real(r8), intent(in)    :: marsh_tidal_range(LBi:,LBj:)
#  if defined MARSH_BIOMASS_VEG
      real(r8), intent(inout) :: plant(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: marsh_biomass_peak(LBi:,LBj:)
      real(r8), intent(inout) :: marsh_vert_rate(LBi:,LBj:)
      real(r8), intent(inout) :: marsh_accret(LBi:,LBj:)
# else
      real(r8), intent(inout) :: marsh_mask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in)    :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(in)    :: marsh_high_water(LBi:UBi,LBj:UBj)
!!#  if defined MARSH_MCKEE_FORMULATION
      real(r8), intent(in)    :: marsh_low_water(LBi:UBi,LBj:UBj)
!!#  endif
!!      real(r8), intent(in)    :: marsh_tidal_range(LBi:UBi,LBj:UBj)
#  if defined MARSH_BIOMASS_VEG
      real(r8), intent(inout) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
#  endif
      real(r8), intent(inout) :: marsh_biomass_peak(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: marsh_vert_rate(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: marsh_accret(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer  :: i, j, k, ised
# if defined MARSH_BIOMASS_VEG
      integer :: ivpr, iveg
# endif
!
      real(r8) :: Dmin, Dmax, depth
      real(r8) :: Samp, Bpeak, Rref, AMC
      real(r8) :: ramp, cff, mtr
      real(r8) :: marsh_bulk_density
      real(r8) :: marsh_vert_rate_insec
      real(r8) :: slope_in_meter,  marsh_organic_matter
!!    real(r8), parameter :: lower_lev = 1.0e-9         ! 800-1000 kg/m3
      real(r8), parameter :: one_over_year2sec = 1.0_r8 /               &
     &                                           (365.0_r8*86400.0_r8)
!
# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Compute marsh vertical growth due to biomass production.
!------------------------------------------------------------------------
!
# if defined WET_DRY && defined MARSH_COLONIZE
      DO j=Jstr,Jend
        DO i=Istr,Iend
!!        marsh_mask(i,j)=marsh_mask(i,j)*rmask_wet(i,j)
          IF (rmask_wet(i,j).eq.0.0_r8) THEN
            marsh_mask(i,j)=1.0_r8
          END IF
        END DO
      END DO
# endif
!
! A ramp could be used to introduce vegetation in the flow.
! Ntimes marsh only activated with marsh tidal range calculation.
!
!   ramp=MIN((FLOOR(tdays(ng)/NTIMES_MARSH)), 1.0_r8 )
!
      J_LOOP: DO j=Jstr,Jend
        DO i=Istr,Iend
!
          IF (tdays(ng).lt.tdays_marsh_growth(ng)) THEN
#  if defined MARSH_KIRWAN_FORMULATION
! Kirwan formulation.
#   if defined MARSH_TIDAL_RANGE_INTERNAL
! Depend on internal calculation.
            Dmin=marsh_high_water(i,j)
! Assume MTR=2*MHW
            mtr=2.0_r8*Dmin
#   else
! Depend on user input.
            mtr=marsh_high_water(i,j)-marsh_low_water(i,j)
#   endif
            Dmax=Par_fac1(ng)*mtr+Par_fac2(ng)+Dmin
#  endif
!
#  if defined MARSH_MCKEE_FORMULATION
! Mckee formulation.
#   if defined MARSH_TIDAL_RANGE_INTERNAL
! Depend on internal calculation.
            Samp=(marsh_high_water(i,j))
#   else
! Depend on user input.
            Samp=0.5*(marsh_high_water(i,j)-marsh_low_water(i,j))
#   endif
            Dmin= 0.429_r8*Samp*2.0_r8+0.253_r8
            Dmax=-0.237_r8*Samp*2.0_r8+0.092_r8
#  endif
! Dmin is higher positive number than Dmax and a higher elevation
! means more towards Dmin and less towards Dmax
!
            Depth=-h(i,j)
!
! scales parabola from 0 to 1 in next line.
!
            cff=0.25_r8*(Dmin-Dmax)*(Dmax-Dmin)
!
! peak biomass as a function of depth.
!
            Bpeak=BMax(ng)*(Depth-Dmax)*(Depth-Dmin)/cff
!
! Not allow negative Bpeak.
!
            Bpeak=MAX(Bpeak,0.0_r8)*marsh_mask(i,j)
!
! Integrated per year amount of below ground biomass.
!  kg/(sq.m-year)
!
            AMC=180.0_r8*Bpeak*nuGp(ng)
          ELSE
!
! 180 growing days, after 180 days make AMC=0, Bpeak=0.0.
!
            Bpeak=0.0_r8
            AMC=0.0_r8
          ENDIF
!
! Masking the amount of marsh biomass and saving in 2D array for I/O.
! (kg/(sq.m-year)
!
          marsh_biomass_peak(i,j)=Bpeak
!
! what material remains after accounting for recalcitrant Carbon.
!  kg/(sq.m-year)
!
          Rref=AMC*chiref(ng)
!
!  marsh_vert_rate is kg/(sq.m-year)/(kg/m^3) = m/year
!    (0.9 is the porosity! )
!
          marsh_bulk_density    =Srho(1,ng)*(1.0_r8-bed(i,j,1,iporo))
          marsh_vert_rate(i,j)  =(Rref/marsh_bulk_density)
!
!  convert m/year to m/s
!
          marsh_vert_rate_insec =marsh_vert_rate(i,j)*one_over_year2sec
!
! convert the rate to m of accretion
!
          slope_in_meter        =marsh_vert_rate_insec*dt(ng)
          marsh_accret(i,j)     =marsh_accret(i,j)+slope_in_meter
!
! Assuming that organic sed in class 1, add bed mass.
! marsh_bulk_density = kg/(m^3) (1-poro)*srho.
!
          marsh_organic_matter  =slope_in_meter*marsh_bulk_density
          bed_mass(i,j,1,nnew,1)=bed_mass(i,j,1,nnew,1)+                &
      &                          marsh_organic_matter
!
#   if defined MARSH_COLONIZE
!
! De-Colonize if marsh biomass production goes to zero
!
!            marsh_mask(i,j)=MIN(1.0_r8,MAX(Bpeak,0.0_r8))
          IF(Bpeak.gt.0.0_r8) THEN
            marsh_mask(i,j)=1.0_r8
          ELSE
            marsh_mask(i,j)=0.0_r8
          END IF
#   endif
#  if defined MARSH_BIOMASS_VEG
!          IF(tdays(ng).gt.tdays_marsh_growth) THEN
!
! Colonize if marsh biomass production (Bpeak)>0.0
!
!          marsh_mask(i,j)=MIN(1.0_r8,MAX(Bpeak,0.0_r8))
!
          DO iveg=1,NVEG
!
! Inside veg loop to change veg properties based on Bpeak
!
            cff=( (Bpeak*1000.0_r8)**beta_pdens(ng) )*marsh_mask(i,j)
            plant(i,j,iveg,isDens)=alpha_pdens(ng)*cff
!
            cff=( (Bpeak*1000.0_r8)**beta_phght(ng) )*marsh_mask(i,j)
            plant(i,j,iveg,isHght)=alpha_phght(ng)*cff
!
            cff=( (Bpeak*1000.0_r8)**beta_pdiam(ng) )*marsh_mask(i,j)
            plant(i,j,iveg,isDiam)=alpha_pdiam(ng)*cff
!
          END DO
!          END IF
#  endif
        END DO
      END DO J_LOOP
!
!---------------------------------------------------------------------
!  Apply laterak boundary conditions (for output purposes only).
!---------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
      END DO
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_biomass_peak)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_vert_rate)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_accret)

# if defined MARSH_BIOMASS_VEG
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_mask)

      DO ivpr=1,NVEGP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, NVEG,                  &
     &                    plant(:,:,:,ivpr))
      END DO
# endif
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_mass(:,:,:,nnew,:))

      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    marsh_biomass_peak,                           &
     &                    marsh_vert_rate,                              &
     &                    marsh_accret)

#  if defined MARSH_BIOMASS_VEG
!
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    marsh_mask)

      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, NVEG, 1, NVEGP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    plant(:,:,:,:))
#  endif
# endif
!
      RETURN
      END SUBROUTINE marsh_vert_growth_tile
!
#endif
      END MODULE marsh_vert_growth_mod

