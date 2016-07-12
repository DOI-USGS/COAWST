#include "cppdefs.h"

      MODULE ad_rho_eos_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_rho_eos.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_rho_eos

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
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
      CALL wclock_on (ng, iADM, 14)
# endif
      CALL ad_rho_eos_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef VAR_RHO_2D
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
# endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % ad_z_r,                          &
# if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % ad_z_w,                          &
# endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t,                           &
# ifdef VAR_RHO_2D
     &                      COUPLING(ng) % ad_rhoA,                     &
     &                      COUPLING(ng) % ad_rhoS,                     &
# endif
# ifdef BV_FREQUENCY_NOT_YET
     &                      MIXING(ng) % ad_bvf,                        &
# endif
# if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
     defined BULK_FLUXES_NOT_YET
     &                      MIXING(ng) % ad_alpha,                      &
     &                      MIXING(ng) % ad_beta,                       &
#  ifdef LMD_DDMIX_NOT_YET
     &                      MIXING(ng) % ad_alfaobeta,                  &
#  endif
# endif
# if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                      OCEAN(ng) % ad_pden,                        &
# endif
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % ad_rho)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 14)
# endif
      RETURN
      END SUBROUTINE ad_rho_eos

# ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE ad_rho_eos_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs,                                 &
#  ifdef MASKING
     &                            rmask,                                &
#  endif
#  ifdef VAR_RHO_2D
     &                            Hz, ad_Hz,                            &
#  endif
     &                            z_r, ad_z_r,                          &
#  if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
     &                            z_w, ad_z_w,                          &
#  endif
     &                            t, ad_t,                              &
#  ifdef VAR_RHO_2D
     &                            ad_rhoA, ad_rhoS,                     &
#  endif
#  ifdef BV_FREQUENCY_NOT_YET
     &                            ad_bvf,                               &
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
     &                            ad_alpha, ad_beta,                    &
#   ifdef LMD_DDMIX_NOT_YET
     &                            ad_alfaobeta,                         &
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                            ad_pden,                              &
#  endif
     &                            rho, ad_rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
#  ifdef SEDIMENT_NOT_YET
      USE mod_sediment
#  endif
!
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#   if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)

#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
#   if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(inout) :: ad_bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
      real(r8), intent(inout) :: ad_alpha(LBi:,LBj:)
      real(r8), intent(inout) :: ad_beta(LBi:,LBj:)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(inout) :: ad_alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(inout) :: ad_pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#   if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))

#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
#   if defined BV_FREQUENCY_NOT_YET || defined VAR_RHO_2D
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(inout) :: ad_bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
      real(r8), intent(inout) :: ad_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(inout) :: ad_alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(inout) :: ad_pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k

      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
      real(r8) :: ad_SedDen, ad_Tp, ad_Tpr10, ad_Ts, ad_Tt, ad_sqrtTs
#  ifdef BV_FREQUENCY_NOT_YET
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
      real(r8) :: ad_bulk_dn, ad_bulk_up, ad_den_dn, ad_den_up
#  endif
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3
      real(r8) :: adfac, adfac1, adfac2, adfac3

      real(r8), dimension(0:9) :: C
      real(r8), dimension(0:9) :: ad_C
#  ifdef EOS_TDERIVATIVE
      real(r8), dimension(0:9) :: dCdT(0:9)
      real(r8), dimension(0:9) :: ad_dCdT(0:9)
      real(r8), dimension(0:9) :: d2Cd2T(0:9)

      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk

      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_wrk
#  endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1

      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_den1

#  ifdef VAR_RHO_2D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhoA1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhoS1
#  endif

#  include "set_bounds.h"
!
!-------------------------------------------------------------------------
!  Initialize adjoint private variables.
!-------------------------------------------------------------------------
!
      ad_Tt=0.0_r8
      ad_Ts=0.0_r8
      ad_Tp=0.0_r8
      ad_Tpr10=0.0_r8
#  ifdef BV_FREQUENCY_NOT_YET
      ad_bulk_dn=0.0_r8
      ad_bulk_up=0.0_r8
      ad_den_dn=0.0_r8
      ad_den_up=0.0_r8
#  endif
      ad_sqrtTs=0.0_r8
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8

      ad_C=0.0_r8
      ad_dCdT=0.0_r8

      DO k=1,N(ng)
        DO i=IminS,ImaxS
#  ifdef EOS_TDERIVATIVE
          ad_DbulkDS(i,k)=0.0_r8
          ad_DbulkDT(i,k)=0.0_r8
          ad_Dden1DS(i,k)=0.0_r8
          ad_Dden1DT(i,k)=0.0_r8
          ad_Scof(i,k)=0.0_r8
          ad_Tcof(i,k)=0.0_r8
          ad_wrk(i,k)=0.0_r8
#  endif
          ad_bulk(i,k)=0.0_r8
          ad_bulk0(i,k)=0.0_r8
          ad_bulk1(i,k)=0.0_r8
          ad_bulk2(i,k)=0.0_r8
          ad_den(i,k)=0.0_r8
          ad_den1(i,k)=0.0_r8
        END DO
      END DO
!
!=======================================================================
!  Adjoint  nonlinear equation of state.  Notice that this equation
!  of state is only valid for potential temperature range of -2C to 40C
!  and a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
#   ifdef BV_FREQUENCY_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_bvf)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_bvf)
#   endif
#   ifdef VAR_RHO_2D
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_rhoA, tl_rhoS)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_rhoA, ad_rhoS)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_alpha, tl_beta)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_alpha, ad_beta)
#    ifdef LMD_DDMIX_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_alfaobeta)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_alfaobeta)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_pden)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_pden)
#   endif
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        NghostPoints,                             &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_rho)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_rho)
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#  ifdef BV_FREQUENCY_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          tl_bvf)
!>
        CALL ad_exchange_w3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                             ad_bvf)
#  endif
#  ifdef VAR_RHO_2D
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_rhoS)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_rhoS)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_rhoA)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_rhoA)
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_beta)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_beta)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_alpha)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_alpha)
#   ifdef LMD_DDMIX_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          tl_alfaobeta)
!>
        CALL ad_exchange_w3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                             ad_alfaobeta)
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_pden)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_pden)
#  endif
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_rho)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_rho)
      END IF
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE related variables.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#  else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#  endif
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!  Compute local nonlinear equation of state coefficients and their
!  derivatives when appropriate.
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
!
            d2Cd2T(0)=2.0_r8*Q02+Tt*(6.0_r8*Q03+Tt*(12.0_r8*Q04+        &
     &                               Tt*20.0_r8*Q05))
            d2Cd2T(1)=2.0_r8*U02+Tt*(6.0_r8*U03+Tt*12.0_r8*U04)
            d2Cd2T(2)=2.0_r8*V02
            d2Cd2T(3)=2.0_r8*A02+Tt*(6.0_r8*A03+Tt*12.0_r8*A04)
            d2Cd2T(4)=2.0_r8*B02+Tt*6.0_r8*B03
            d2Cd2T(5)=2.0_r8*D02
            d2Cd2T(6)=2.0_r8*E02+Tt*6.0_r8*E03
            d2Cd2T(7)=2.0_r8*F02
            d2Cd2T(8)=2.0_r8*G03
            d2Cd2T(9)=2.0_r8*H02
#  endif
!
!  Compute BASIC STATE density (kg/m3) at standard one atmosphere
!  pressure.
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

#  ifdef EOS_TDERIVATIVE
!
!  Compute BASIC STATE d(den1)/d(S) and d(den1)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            Dden1DS(i,k)=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT(i,k)=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
#  endif
!
!  Compute BASIC STATE secant bulk modulus.
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!
!  Compute BASIC STATE d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            DbulkDS(i,k)=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT(i,k)=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
#  endif
!
!  Compute local "in situ" density anomaly (kg/m3 - 1000).  The (i,k)
!  DO-loop is closed here because of the adjoint to facilitate vertical
!  integrals of the BASIC STATE.
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
#  if defined SEDIMENT_NOT_YET && defined SED_DENS_NOT_YET
            SedDen=0.0_r8
            DO ised=1,NST
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,idsed(ised))*                         &
     &               (Srho(ised,ng)-den(i,k))*cff1
            END DO
            den(i,k)=den(i,k)+SedDen
#  endif
            den(i,k)=den(i,k)-1000.0_r8
#  ifdef MASKING
            den(i,k)=den(i,k)*rmask(i,j)
#  endif
          END DO
        END DO

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE thermal expansion (1/Celsius) and saline
!  contraction (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX_NOT_YET
        DO k=1,N(ng)
#   else
        DO k=N(ng),N(ng)
#   endif
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
          END DO
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint "in situ" density anomaly (kg/m3 - 1000) and adjoint
!  potential density anomaly (kg/m3 - 1000) referenced to the surface
!  into global arrays.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
#   ifdef MASKING
!>          tl_pden(i,j,k)=tl_pden(i,k)*rmask(i,j)
!>
            ad_pden(i,j,k)=ad_pden(i,k)*rmask(i,j)
#   endif
!>          tl_pden(i,j,k)=tl_den1(i,k)             ! This gives a fatal
!>                                                  ! result in 4D-Var
            ad_den1(i,k)=ad_den1(i,k)+ad_pden(i,j,k)! posterior error...
            ad_pden(i,j,k)=0.0_r8
#  endif
!>          tl_rho(i,j,k)=tl_den(i,k)
!>
            ad_den(i,k)=ad_den(i,k)+ad_rho(i,j,k)
            ad_rho(i,j,k)=0.0_r8
          END DO
        END DO

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute adjoint thermal expansion (1/Celsius) and adjoint saline
!  contraction (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX_NOT_YET
!>      DO k=1,N(ng)
!>
        DO k=N(ng),1,-1
#   else
        DO k=N(ng),N(ng)
#   endif
          IF (k.eq.N(ng)) THEN
            DO i=IstrT,IendT
              cff=1.0_r8/wrk(i,N(ng))
!>            tl_beta (i,j)=tl_cff*Scof(i,N(ng))+cff*tl_Scof(i,N(ng))
!>            tl_alpha(i,j)=tl_cff*Tcof(i,N(ng))+cff*tl_Tcof(i,N(ng))
!>
              ad_Scof(i,N(ng))=ad_Scof(i,N(ng))+cff*ad_beta (i,j)
              ad_Tcof(i,N(ng))=ad_Tcof(i,N(ng))+cff*ad_alpha(i,j)
              ad_cff=ad_beta (i,j)*Scof(i,N(ng))+                       &
     &               ad_alpha(i,j)*Tcof(i,N(ng))
              ad_beta (i,j)=0.0_r8
              ad_alpha(i,j)=0.0_r8
!>            tl_cff=-cff*cff*tl_wrk(i,N(ng))
!>
              ad_wrk(i,N(ng))=ad_wrk(i,N(ng))-cff*cff*ad_cff
              ad_cff=0.0_r8
            END DO
          END IF
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
#   ifdef LMD_DDMIX_NOT_YET
!>          tl_alfaobeta(i,j,k)=(tl_Tcof(i,k)*Scof(i,k)-                &
!>   &                           Tcof(i,k)*tl_Scof(i,k))/               &
!>   &                          (Scof(i,k)*Scof(i,k))
!>
            adfac=ad_alfaobeta(i,j,k)/(Scof(i,k)*Scof(i,k))
            ad_Tcof(i,k)=ad_Tcof(i,k)+Scof(i,k)*adfac
            ad_Scof(i,k)=ad_Scof(i,k)-Tcof(i,k)*adfac
            ad_alfaobeta(i,j,k)=0.0_r8
#   endif
!>          tl_Scof(i,k)= (tl_DbulkDS(i,k)*cff1+                        &
!>   &                     DbulkDS(i,k)*tl_cff1+                        &
!>   &                     tl_Dden1DS(i,k)*cff2+                        &
!>   &                     Dden1DS(i,k)*tl_cff2)
!>
            ad_DbulkDS(i,k)=ad_DbulkDS(i,k)+ad_Scof(i,k)*cff1
            ad_Dden1DS(i,k)=ad_Dden1DS(i,k)+ad_Scof(i,k)*cff2
            ad_cff1=DbulkDS(i,k)*ad_Scof(i,k)
            ad_cff2=Dden1DS(i,k)*ad_Scof(i,k)
            ad_Scof(i,k)=0.0_r8
!>          tl_Tcof(i,k)=-(tl_DbulkDT(i,k)*cff1+                        &
!>   &                     DbulkDT(i,k)*tl_cff1+                        &
!>   &                     tl_Dden1DT(i,k)*cff2+                        &
!>   &                     Dden1DT(i,k)*tl_cff2)
!>
            ad_DbulkDT(i,k)=ad_DbulkDT(i,k)-ad_Tcof(i,k)*cff1
            ad_Dden1DT(i,k)=ad_Dden1DT(i,k)-ad_Tcof(i,k)*cff2
            ad_cff1=ad_cff1-DbulkDT(i,k)*ad_Tcof(i,k)
            ad_cff2=ad_cff2-Dden1DT(i,k)*ad_Tcof(i,k)
            ad_Tcof(i,k)=0.0_r8
!>          tl_wrk(i,k)=cff*(cff*tl_den(i,k)+                           &
!>   &                       2.0_r8*tl_cff*(den(i,k)+1000.0_r8))
!>
            adfac=cff*ad_wrk(i,k)
            ad_den(i,k)=ad_den(i,k)+cff*adfac
            ad_cff=ad_cff+2.0_r8*(den(i,k)+1000.0_r8)*adfac
            ad_wrk(i,k)=0.0_r8
!>          tl_cff2=tl_bulk(i,k)*cff+bulk(i,k)*tl_cff
!>
            ad_bulk(i,k)=ad_bulk(i,k)+ad_cff2*cff
            ad_cff=ad_cff+bulk(i,k)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_Tpr10*den1(i,k)+Tpr10*tl_den1(i,k)
!>
            ad_Tpr10=ad_Tpr10+ad_cff1*den1(i,k)
            ad_den1(i,k)=ad_den1(i,k)+Tpr10*ad_cff1
            ad_cff1=0.0_r8
!>          tl_cff=tl_bulk(i,k)+tl_Tpr10
!>
            ad_bulk(i,k)=ad_bulk(i,k)+ad_cff
            ad_Tpr10=ad_Tpr10+ad_cff
            ad_cff=0.0_r8
!>          tl_Tpr10=0.1_r8*tl_z_r(i,j,k)
!>
            ad_z_r(i,j,k)=ad_z_r(i,j,k)+0.1_r8*ad_Tpr10
            ad_Tpr10=0.0_r8
          END DO
        END DO
#  endif

#  if defined BV_FREQUENCY_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute adjoint Brunt-Vaisala frequency (1/s2) at horizontal
!  RHO-points and vertical W-points.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
!>        tl_bvf(i,j,N(ng))=0.0_r8
!>
          ad_bvf(i,j,N(ng))=0.0_r8
!>        tl_bvf(i,j,0)=0.0_r8
!>
          ad_bvf(i,j,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
            cff3=1.0_r8/(0.5_r8*(den_up+den_dn)*                        &
     &                   (z_r(i,j,k+1)-z_r(i,j,k)))
!>          tl_bvf(i,j,k)=-g*((tl_den_up-tl_den_dn)*cff3+               &
!>   &                        (den_up-den_dn)*tl_cff3)
!>
            adfac=-g*ad_bvf(i,j,k)
            adfac1=adfac*cff3
            ad_cff3=ad_cff3+(den_up-den_dn)*adfac
            ad_den_up=ad_den_up+adfac1
            ad_den_dn=ad_den_dn-adfac1
            ad_bvf(i,j,k)=0.0_r8
!>          tl_cff3=-cff3*cff3*                                         &
!>   &               0.5_r8*((tl_den_up+tl_den_dn)*                     &
!>   &                       (z_r(i,j,k+1)-z_r(i,j,k))+                 &
!>   &                       (den_up+den_dn)*                           &
!>   &                       (tl_z_r(i,j,k+1)-tl_z_r(i,j,k)))
!>
            adfac=-cff3*cff3*0.5_r8*ad_cff3
            adfac1=adfac*(z_r(i,j,k+1)-z_r(i,j,k))
            adfac2=adfac*(den_up+den_dn)
            ad_z_r(i,j,k  )=ad_z_r(i,j,k  )-adfac2
            ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)+adfac2
            ad_den_up=ad_den_up+adfac1
            ad_den_dn=ad_den_dn+adfac1
            ad_cff3=0.0_r8
!>          tl_den_dn=tl_cff2*(den1(i,k  )*bulk_dn)+                    &
!>   &                cff2*(tl_den1(i,k  )*bulk_dn+                     &
!>   &                      den1(i,k  )*tl_bulk_dn)
!>          tl_den_up=tl_cff1*(den1(i,k+1)*bulk_up)+                    &
!>   &                cff1*(tl_den1(i,k+1)*bulk_up+                     &
!>   &                      den1(i,k+1)*tl_bulk_up)
!>
            adfac1=cff2*ad_den_dn
            adfac2=cff1*ad_den_up
            ad_cff2=ad_cff2+(den1(i,k  )*bulk_dn)*ad_den_dn
            ad_cff1=ad_cff1+(den1(i,k+1)*bulk_up)*ad_den_up
            ad_den1(i,k  )=ad_den1(i,k  )+bulk_dn*adfac1
            ad_den1(i,k+1)=ad_den1(i,k+1)+bulk_up*adfac2
            ad_bulk_dn=ad_bulk_dn+den1(i,k  )*adfac1
            ad_bulk_up=ad_bulk_up+den1(i,k+1)*adfac2
            ad_den_dn=0.0_r8
            ad_den_up=0.0_r8
!>          tl_cff2=-cff2*cff2*(tl_bulk_dn+0.1_r8*tl_z_w(i,j,k))
!>          tl_cff1=-cff1*cff1*(tl_bulk_up+0.1_r8*tl_z_w(i,j,k))
!>
            adfac1=-cff2*cff2*ad_cff2
            adfac2=-cff1*cff1*ad_cff1
            ad_bulk_dn=ad_bulk_dn+adfac1
            ad_bulk_up=ad_bulk_up+adfac2
            ad_z_w(i,j,k)=ad_z_w(i,j,k)+                                &
     &                    0.1_r8*(adfac1+adfac2)
            ad_cff2=0.0_r8
            ad_cff1=0.0_r8
!>          tl_bulk_dn=tl_bulk0(i,k  )-                                 &
!>   &                 tl_z_w(i,j,k)*(bulk1(i,k  )-                     &
!>   &                                bulk2(i,k  )*z_w(i,j,k))-         &
!>   &                 z_w(i,j,k)*(tl_bulk1(i,k  )-                     &
!>   &                             tl_bulk2(i,k  )*z_w(i,j,k)-          &
!>   &                             bulk2(i,k  )*tl_z_w(i,j,k))
!>          tl_bulk_up=tl_bulk0(i,k+1)-                                 &
!>   &                 tl_z_w(i,j,k)*(bulk1(i,k+1)-                     &
!>   &                                bulk2(i,k+1)*z_w(i,j,k))-         &
!>   &                 z_w(i,j,k)*(tl_bulk1(i,k+1)-                     &
!>   &                             tl_bulk2(i,k+1)*z_w(i,j,k)-          &
!>   &                             bulk2(i,k+1)*tl_z_w(i,j,k))
!>
            adfac1=z_w(i,j,k)*ad_bulk_dn
            adfac2=z_w(i,j,k)*ad_bulk_up
            ad_bulk0(i,k  )=ad_bulk0(i,k  )+ad_bulk_dn
            ad_bulk0(i,k+1)=ad_bulk0(i,k+1)+ad_bulk_up
            ad_z_w(i,j,k)=ad_z_w(i,j,k)-                                &
     &                    (bulk1(i,k  )-                                &
     &                     bulk2(i,k  )*z_w(i,j,k)-                     &
     &                     bulk2(i,k  ))*ad_bulk_dn-                    &
     &                    (bulk1(i,k+1)-                                &
     &                     bulk2(i,k+1)*z_w(i,j,k)-                     &
     &                     bulk2(i,k+1))*ad_bulk_up
            ad_bulk1(i,k  )=ad_bulk1(i,k  )-adfac1
            ad_bulk1(i,k+1)=ad_bulk1(i,k+1)-adfac2
            ad_bulk2(i,k  )=ad_bulk2(i,k  )+z_w(i,j,k)*adfac1
            ad_bulk2(i,k+1)=ad_bulk2(i,k+1)+z_w(i,j,k)*adfac2
            ad_bulk_dn=0.0_r8
            ad_bulk_up=0.0_r8
          END DO
        END DO
#  endif

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute adjoint vertical averaged density (ad_rhoA) and adjoint
!  density perturbation (ad_rhoS) used in adjoint barotropic pressure
!  gradient.
!-----------------------------------------------------------------------
!
!  Compute temporary intermediate BASIC STATE "rhoS1" and "rhoA1".
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS1(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA1(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            rhoS1(i,j)=rhoS1(i,j)+Hz(i,j,k)*(rhoA1(i,j)+0.5_r8*cff1)
            rhoA1(i,j)=rhoA1(i,j)+cff1
          END DO
        END DO
!
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
!>        tl_rhoS(i,j)=2.0_r8*cff2*                                     &
!>   &                 cff1*(2.0_r8*tl_cff1*rhoS1(i,j)+                 &
!>   &                       cff1*tl_rhoS(i,j))
!>
          adfac=2.0_r8*cff2*cff1*ad_rhoS(i,j)
          ad_cff1=2.0_r8*rhoS1(i,j)*adfac
          ad_rhoS(i,j)=cff1*adfac
!>        tl_rhoA(i,j)=cff2*(tl_cff1*rhoA1(i,j)+cff1*tl_rhoA(i,j))
!>
          adfac=cff2*ad_rhoA(i,j)
          ad_cff1=ad_cff1+rhoA1(i,j)*adfac
          ad_rhoA(i,j)=cff1*adfac
!>        tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
!>
          adfac=-cff1*cff1*ad_cff1
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac
          ad_z_w(i,j,0    )=ad_z_w(i,j,0    )-adfac
          ad_cff1=0.0_r8
        END DO
!
!  Compute appropriate intermediate BASIC STATE "rhoA1".
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoA1(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
!>          tl_rhoA(i,j)=tl_rhoA(i,j)+tl_cff1
!>
            ad_cff1=ad_rhoA(i,j)
!>          tl_rhoS(i,j)=tl_rhoS(i,j)+                                  &
!>   &                   tl_Hz(i,j,k)*(rhoA1(i,j)+0.5_r8*cff1)+         &
!>   &                   Hz(i,j,k)*(tl_rhoA(i,j)+0.5_r8*tl_cff1)
!>
            adfac=Hz(i,j,k)*ad_rhoS(i,j)
            ad_rhoA(i,j)=ad_rhoA(i,j)+adfac
            ad_cff1=ad_cff1+0.5_r8*adfac
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                  &
     &                   (rhoA1(i,j)+0.5_r8*cff1)*ad_rhoS(i,j)
!>          tl_cff1=tl_den(i,k)*Hz(i,j,k)+                              &
!>   &              den(i,k)*tl_Hz(i,j,k)
!>
            ad_den(i,k)=ad_den(i,k)+Hz(i,j,k)*ad_cff1
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+den(i,k)*ad_cff1
            ad_cff1=0.0_r8
            rhoA1(i,j)=rhoA1(i,j)+cff1
          END DO
        END DO
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
!>        tl_rhoA(i,j)=tl_cff1
!>
          ad_cff1=ad_rhoA(i,j)
          ad_rhoA(i,j)=0.0_r8
!>        tl_rhoS(i,j)=0.5_r8*(tl_cff1*Hz(i,j,N(ng))+                   &
!>   &                         cff1*tl_Hz(i,j,N(ng)))
!>
          adfac=0.5_r8*ad_rhoS(i,j)
          ad_cff1=ad_cff1+Hz(i,j,N(ng))*adfac
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+cff1*adfac
          ad_rhoS(i,j)=0.0_r8
!>        tl_cff1=tl_den(i,N(ng))*Hz(i,j,N(ng))+                        &
!>   &            den(i,N(ng))*tl_Hz(i,j,N(ng))
!>
          ad_den(i,N(ng))=ad_den(i,N(ng))+Hz(i,j,N(ng))*ad_cff1
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+den(i,N(ng))*ad_cff1
          ad_cff1=0.0_r8
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Adjoint nonlinear equation of state.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
!
!  Check temperature and salinity minimum valid values. Assign depth
!  to the pressure.
!
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#  else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#  endif
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!  Compute local nonlinear equation of state coefficients and their
!  derivatives when appropriate.  These coefficients can be stored
!  in slab (i,k) arrays to avoid recompute them twice. However, the
!  equivalent of 50 slabs arrays are required.
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
!
            d2Cd2T(0)=2.0_r8*Q02+Tt*(6.0_r8*Q03+Tt*(12.0_r8*Q04+        &
     &                               Tt*20.0_r8*Q05))
            d2Cd2T(1)=2.0_r8*U02+Tt*(6.0_r8*U03+Tt*12.0_r8*U04)
            d2Cd2T(2)=2.0_r8*V02
            d2Cd2T(3)=2.0_r8*A02+Tt*(6.0_r8*A03+Tt*12.0_r8*A04)
            d2Cd2T(4)=2.0_r8*B02+Tt*6.0_r8*B03
            d2Cd2T(5)=2.0_r8*D02
            d2Cd2T(6)=2.0_r8*E02+Tt*6.0_r8*E03
            d2Cd2T(7)=2.0_r8*F02
            d2Cd2T(8)=2.0_r8*G03
            d2Cd2T(9)=2.0_r8*H02
#  endif
!
!-----------------------------------------------------------------------
!  Compute local adjoint "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
#  ifdef MASKING
!>          tl_den(i,k)=tl_den(i,k)*rmask(i,j)
!>
            ad_den(i,k)=ad_den(i,k)*rmask(i,j)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_DENS_NOT_YET
!>          tl_den(i,k)=tl_den(i,k)+tl_SedDen
!>
            ad_SedDen=ad_SedDen+ad_den(i,k)
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
!>            tl_SedDen=tl_SedDen+                                      &
!>   &                  (tl_t(i,j,k,nrhs,itrc)*                         &
!>   &                   (Srho(ised,ng)-den(i,k))-                      &
!>   &                   t(i,j,k,nrhs,itrc)*                            &
!>   &                   tl_den(i,k))*cff1
!>
              adfac=cff1*ad_SedDen
              ad_den(i,k)=ad_den(i,k)-                                  &
     &                    t(i,j,k,nrhs,idsed(ised))*adfac
              ad_t(i,j,k,nrhs,itrc)=ad_t(i,j,k,nrhs,itrc)+              &
     &                              (Srho(ised,ng)-den(i,k))*adfac
            END DO
!>          tl_SedDen=0.0_r8
!>
            ad_SedDen=0.0_r8
#  endif
!>          tl_den(i,k)=tl_den1(i,k)*bulk(i,k)*cff+                     &
!>   &                  den1(i,k)*(tl_bulk(i,k)*cff+                    &
!>   &                             bulk(i,k)*tl_cff)
!>
            adfac1=den1(i,k)*ad_den(i,k)
            ad_den1(i,k)=ad_den1(i,k)+bulk(i,k)*cff*ad_den(i,k)
            ad_bulk(i,k)=ad_bulk(i,k)+cff*adfac1
            ad_cff=ad_cff+bulk(i,k)*adfac1
            ad_den(i,k)=0.0_r8
!>          tl_cff=-cff*cff*(tl_bulk(i,k)+tl_Tpr10)
!>
            adfac=-cff*cff*ad_cff
            ad_bulk(i,k)=ad_bulk(i,k)+adfac
            ad_Tpr10=ad_Tpr10+adfac
            ad_cff=0.0_r8

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
!>          tl_DbulkDT(i,k)=tl_dCdT(3)+                                 &
!>   &                      tl_Ts*(dCdT(4)+sqrtTs*dCdT(5))+             &
!>   &                      Ts*(tl_dCdT(4)+tl_sqrtTs*dCdT(5)+           &
!>   &                                     sqrtTs*tl_dCdT(5))-          &
!>   &                      tl_Tp*(dCdT(6)+Ts*dCdT(7)-                  &
!>   &                             Tp*(dCdT(8)+Ts*dCdT(9)))-            &
!>   &                      Tp*(tl_dCdT(6)+tl_Ts*dCdT(7)+Ts*tl_dCdT(7)- &
!>   &                          tl_Tp*(dCdT(8)+Ts*dCdT(9))-             &
!>   &                          Tp*(tl_dCdT(8)+tl_Ts*dCdT(9)+           &
!>   &                                         Ts*tl_dCdT(9)))
!>
            adfac1=Ts*ad_DbulkDT(i,k)
            adfac2=Tp*ad_DbulkDT(i,k)
            adfac3=adfac2*Tp
            ad_dCdT(3)=ad_dCdT(3)+ad_DbulkDT(i,k)
            ad_dCdT(4)=ad_dCdT(4)+adfac1
            ad_dCdT(5)=ad_dCdT(5)+sqrtTs*adfac1
            ad_dCdT(6)=ad_dCdT(6)-adfac2
            ad_dCdT(7)=ad_dCdT(7)-Ts*adfac2
            ad_dCdT(8)=ad_dCdT(8)+adfac3
            ad_dCdT(9)=ad_dCdT(9)+Ts*adfac3
            ad_sqrtTs=ad_sqrtTs+dCdT(5)*adfac1
            ad_Ts=ad_Ts+                                                &
     &            ad_DbulkDT(i,k)*                                      &
     &            (dCdT(4)+sqrtTs*dCdT(5)-                              &
     &             Tp*(dCdT(7)-Tp*dCdT(9)))
            ad_Tp=ad_Tp-                                                &
     &            ad_DbulkDT(i,k)*                                      &
     &            (dCdT(6)+Ts*dCdT(7)-                                  &
     &             2.0_r8*Tp*(dCdT(8)+Ts*dCdT(9)))
            ad_DbulkDT(i,k)=0.0_r8
!>          tl_DbulkDS(i,k)=tl_C(4)+                                    &
!>   &                      1.5_r8*(tl_sqrtTs*C(5)+                     &
!>   &                              sqrtTs*tl_C(5))-                    &
!>   &                      tl_Tp*(C(7)+sqrtTs*1.5_r8*G00-              &
!>   &                             Tp*C(9))-                            &
!>   &                      Tp*(tl_C(7)+tl_sqrtTs*1.5_r8*G00-           &
!>   &                          tl_Tp*C(9)-Tp*tl_C(9))
!>
            adfac1=1.5_r8*ad_DbulkDS(i,k)
            adfac2=Tp*ad_DbulkDS(i,k)
            ad_C(4)=ad_C(4)+ad_DbulkDS(i,k)
            ad_C(5)=ad_C(5)+sqrtTs*adfac1
            ad_C(7)=ad_C(7)-adfac2
            ad_C(9)=ad_C(9)+Tp*adfac2
            ad_sqrtTs=ad_sqrtTs+                                        &
     &                (C(5)-Tp*G00)*adfac1
            ad_Tp=ad_Tp-                                                &
     &            ad_DbulkDS(i,k)*                                      &
     &            (C(7)+sqrtTs*1.5_r8*G00-Tp*C(9)-                      &
     &             Tp*C(9))
            ad_DbulkDS(i,k)=0.0_r8
!>          tl_dCdT(9)=tl_Tt*d2Cd2T(9)
!>          tl_dCdT(8)=tl_Tt*d2Cd2T(8)
!>          tl_dCdT(7)=tl_Tt*d2Cd2T(7)
!>          tl_dCdT(6)=tl_Tt*d2Cd2T(6)
!>          tl_dCdT(5)=tl_Tt*d2Cd2T(5)
!>          tl_dCdT(4)=tl_Tt*d2Cd2T(4)
!>          tl_dCdT(3)=tl_Tt*d2Cd2T(3)
!>
            ad_Tt=ad_Tt+d2Cd2T(9)*ad_dCdT(9)+                           &
     &                  d2Cd2T(8)*ad_dCdT(8)+                           &
     &                  d2Cd2T(7)*ad_dCdT(7)+                           &
     &                  d2Cd2T(6)*ad_dCdT(6)+                           &
     &                  d2Cd2T(5)*ad_dCdT(5)+                           &
     &                  d2Cd2T(4)*ad_dCdT(4)+                           &
     &                  d2Cd2T(3)*ad_dCdT(3)
            ad_dCdT(9)=0.0_r8
            ad_dCdT(8)=0.0_r8
            ad_dCdT(7)=0.0_r8
            ad_dCdT(6)=0.0_r8
            ad_dCdT(5)=0.0_r8
            ad_dCdT(4)=0.0_r8
            ad_dCdT(3)=0.0_r8
#  endif
!
!  Compute adjoint secant bulk modulus.
!
!>          tl_bulk (i,k)=tl_bulk0(i,k)-                                &
!>   &                     tl_Tp*(bulk1(i,k)-Tp*bulk2(i,k))-            &
!>   &                     Tp*(tl_bulk1(i,k)-                           &
!>   &                         tl_Tp*bulk2(i,k)-                        &
!>   &                         Tp*tl_bulk2(i,k))
!>
            adfac=Tp*ad_bulk(i,k)
            ad_bulk0(i,k)=ad_bulk0(i,k)+ad_bulk(i,k)
            ad_bulk1(i,k)=ad_bulk1(i,k)-adfac
            ad_bulk2(i,k)=ad_bulk2(i,k)+adfac*Tp
            ad_Tp=ad_Tp-                                                &
     &            ad_bulk(i,k)*(bulk1(i,k)-Tp*bulk2(i,k))+              &
     &            adfac*bulk2(i,k)
            ad_bulk(i,k)=0.0_r8
!>          tl_bulk2(i,k)=tl_C(8)+tl_Ts*C(9)+Ts*tl_C(9)
!>
            ad_C(8)=ad_C(8)+ad_bulk2(i,k)
            ad_C(9)=ad_C(9)+Ts*ad_bulk2(i,k)
            ad_Ts=ad_Ts+ad_bulk2(i,k)*C(9)
            ad_bulk2(i,k)=0.0_r8
!>          tl_bulk1(i,k)=tl_C(6)+                                      &
!>   &                    tl_Ts*(C(7)+sqrtTs*G00)+                      &
!>   &                    Ts*(tl_C(7)+tl_sqrtTs*G00)
!>
            adfac=Ts*ad_bulk1(i,k)
            ad_C(6)=ad_C(6)+ad_bulk1(i,k)
            ad_C(7)=ad_C(7)+adfac
            ad_Ts=ad_Ts+ad_bulk1(i,k)*(C(7)+sqrtTs*G00)
            ad_sqrtTs=ad_sqrtTs+adfac*G00
            ad_bulk1(i,k)=0.0_r8
!>          tl_bulk0(i,k)=tl_C(3)+                                      &
!>   &                    tl_Ts*(C(4)+sqrtTs*C(5))+                     &
!>   &                    Ts*(tl_C(4)+tl_sqrtTs*C(5)+                   &
!>   &                        sqrtTs*tl_C(5))
!>
            adfac=Ts*ad_bulk0(i,k)
            ad_C(3)=ad_C(3)+ad_bulk0(i,k)
            ad_C(4)=ad_C(4)+adfac
            ad_C(5)=ad_C(5)+sqrtTs*adfac
            ad_Ts=ad_Ts+ad_bulk0(i,k)*(C(4)+sqrtTs*C(5))
            ad_sqrtTs=ad_sqrtTs+C(5)*adfac
            ad_bulk0(i,k)=0.0_r8
!>          tl_C(9)=tl_Tt*dCdT(9)
!>          tl_C(8)=tl_Tt*dCdT(8)
!>          tl_C(7)=tl_Tt*dCdT(7)
!>          tl_C(6)=tl_Tt*dCdT(6)
!>          tl_C(5)=tl_Tt*dCdT(5)
!>          tl_C(4)=tl_Tt*dCdT(4)
!>          tl_C(3)=tl_Tt*dCdT(3)
!>
            ad_Tt=ad_Tt+ad_C(9)*dCdT(9)+                                &
     &                  ad_C(8)*dCdT(8)+                                &
     &                  ad_C(7)*dCdT(7)+                                &
     &                  ad_C(6)*dCdT(6)+                                &
     &                  ad_C(5)*dCdT(5)+                                &
     &                  ad_C(4)*dCdT(4)+                                &
     &                  ad_C(3)*dCdT(3)
            ad_C(9)=0.0_r8
            ad_C(8)=0.0_r8
            ad_C(7)=0.0_r8
            ad_C(6)=0.0_r8
            ad_C(5)=0.0_r8
            ad_C(4)=0.0_r8
            ad_C(3)=0.0_r8

#  ifdef EOS_TDERIVATIVE
!
!  Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
!  computation of thermal expansion and saline contraction
!  coefficients.
!
!>          tl_Dden1DT(i,k)=tl_dCdT(0)+                                 &
!>   &                      tl_Ts*(dCdT(1)+sqrtTs*dCdT(2))+             &
!>   &                      Ts*(tl_dCdT(1)+tl_sqrtTs*dCdT(2)+           &
!>   &                                     sqrtTs*tl_dCdT(2))
!>
            adfac1=Ts*ad_Dden1DT(i,k)
            ad_dCdT(0)=ad_dCdT(0)+ad_Dden1DT(i,k)
            ad_dCdT(1)=ad_dCdT(1)+adfac1
            ad_dCdT(2)=ad_dCdT(2)+sqrtTs*adfac1
            ad_Ts=ad_Ts+                                                &
     &            (dCdT(1)+sqrtTs*dCdT(2))*ad_Dden1DT(i,k)
            ad_sqrtTs=ad_sqrtTs+dCdT(2)*adfac1
            ad_Dden1DT(i,k)=0.0_r8
!>          tl_Dden1DS(i,k)=tl_C(1)+                                    &
!>   &                      1.5_r8*(tl_C(2)*sqrtTs+                     &
!>   &                              C(2)*tl_sqrtTs)+                    &
!>   &                      2.0_r8*W00*tl_Ts
!>
            adfac1=1.5_r8*ad_Dden1DS(i,k)
            ad_C(1)=ad_C(1)+ad_Dden1DS(i,k)
            ad_C(2)=ad_C(2)+sqrtTs*adfac1
            ad_Ts=ad_Ts+2.0_r8*W00*ad_Dden1DS(i,k)
            ad_sqrtTs=ad_sqrtTs+C(2)*adfac1
            ad_Dden1DS(i,k)=0.0_r8
!>          tl_dCdT(2)=tl_Tt*d2Cd2T(2)
!>          tl_dCdT(1)=tl_Tt*d2Cd2T(1)
!>          tl_dCdT(0)=tl_Tt*d2Cd2T(0)
!>
            ad_Tt=ad_Tt+d2Cd2T(2)*ad_dCdT(2)+                           &
     &                  d2Cd2T(1)*ad_dCdT(1)+                           &
     &                  d2Cd2T(0)*ad_dCdT(0)
            ad_dCdT(2)=0.0_r8
            ad_dCdT(1)=0.0_r8
            ad_dCdT(0)=0.0_r8
#  endif
!
!  Compute basic state and tangent linear density (kg/m3) at standard
!  one atmosphere pressure.
!
!>          tl_den1(i,k)=tl_C(0)+                                       &
!>   &                   tl_Ts*(C(1)+sqrtTs*C(2)+Ts*W00)+               &
!>   &                   Ts*(tl_C(1)+tl_sqrtTs*C(2)+                    &
!>   &                       sqrtTs*tl_C(2)+tl_Ts*W00)
!>
            adfac=Ts*ad_den1(i,k)
            ad_C(0)=ad_C(0)+ad_den1(i,k)
            ad_C(1)=ad_C(1)+adfac
            ad_C(2)=ad_C(2)+adfac*sqrtTs
            ad_Ts=ad_Ts+                                                &
     &            ad_den1(i,k)*(C(1)+sqrtTs*C(2)+Ts*W00)+               &
     &            adfac*W00
            ad_sqrtTs=ad_sqrtTs+adfac*C(2)
            ad_den1(i,k)=0.0_r8
!>          tl_C(2)=tl_Tt*dCdT(2)
!>          tl_C(1)=tl_Tt*dCdT(1)
!>          tl_C(0)=tl_Tt*dCdT(0)
!>
            ad_Tt=ad_Tt+ad_C(2)*dCdT(2)+                                &
     &                  ad_C(1)*dCdT(1)+                                &
     &                  ad_C(0)*dCdT(0)
            ad_C(2)=0.0_r8
            ad_C(1)=0.0_r8
            ad_C(0)=0.0_r8
!
!  Check temperature and salinity minimum valid values. Assign depth
!  to the pressure.
!
!>          tl_Tpr10=0.1_r8*tl_Tp
!>
            ad_Tp=ad_Tp+0.1_r8*ad_Tpr10
            ad_Tpr10=0.0_r8
!>          tl_Tp=tl_z_r(i,j,k)
!>
            ad_z_r(i,j,k)=ad_z_r(i,j,k)+ad_Tp
            ad_Tp=0.0_r8

#  ifdef SALINITY
            IF (Ts.ne.0.0_r8) THEN
!>            tl_sqrtTs=0.5_r8*tl_Ts/SQRT(Ts)
!>
              ad_Ts=ad_Ts+0.5_r8*ad_sqrtTs/SQRT(Ts)
              ad_sqrtTs=0.0_r8
            ELSE
!>            tl_sqrtTs=0.0_r8
!>
              ad_sqrtTs=0.0_r8
            END IF
!>          tl_Ts=(0.5_r8-SIGN(0.5_r8,-t(i,j,k,nrhs,isalt)))*
!>   &            tl_t(i,j,k,nrhs,isalt)
!>
            ad_t(i,j,k,nrhs,isalt)=ad_t(i,j,k,nrhs,isalt)+              &
     &                             (0.5_r8-SIGN(0.5_r8,                 &
     &                                          -t(i,j,k,nrhs,isalt)))* &
     &                             ad_Ts
            ad_Ts=0.0_r8
#  else
!>          tl_sqrtTs=0.0_r8
!>
            ad_sqrtTs=0.0_r8
!>          tl_Ts=0.0_r8
!>
            ad_Ts=0.0_r8
#  endif
!>          tl_Tt=(0.5_r8-SIGN(0.5_r8,-2.0_r8-t(i,j,k,nrhs,itemp)))*
!>   &            tl_t(i,j,k,nrhs,itemp)
!>
            ad_t(i,j,k,nrhs,itemp)=ad_t(i,j,k,nrhs,itemp)+              &
     &                             (0.5_r8-SIGN(0.5_r8,-2.0_r8-         &
     &                                          t(i,j,k,nrhs,itemp)))*  &
     &                             ad_Tt
            ad_Tt=0.0_r8
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE ad_rho_eos_tile
# endif

# ifndef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE ad_rho_eos_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs,                                 &
#  ifdef MASKING
     &                            rmask,                                &
#  endif
#  ifdef VAR_RHO_2D
     &                            Hz, ad_Hz,                            &
#  endif
     &                            z_r, ad_z_r,                          &
#  ifdef VAR_RHO_2D
     &                            z_w, ad_z_w,                          &
#  endif
     &                            t, ad_t,                              &
#  ifdef VAR_RHO_2D
     &                            ad_rhoA, ad_rhoS,                     &
#  endif
#  ifdef BV_FREQUENCY_NOT_YET
     &                            ad_bvf,                               &
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
     &                            ad_alpha, ad_beta,                    &
#   ifdef LMD_DDMIX_NOT_YET
     &                            ad_alfaobeta,                         &
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                            ad_pden,                              &
#  endif
     &                            rho, ad_rho)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
#  ifdef SEDIMENT_NOT_YET
      USE mod_sediment
#  endif
!
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(inout) :: tl_bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
      real(r8), intent(inout) :: ad_alpha(LBi:,LBj:)
      real(r8), intent(inout) :: ad_beta(LBi:,LBj:)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(inout) :: ad_alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(inout) :: ad_pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

      real(r8), intent(inout) :: rho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(inout) :: ad_bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
      real(r8), intent(inout) :: ad_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(inout) :: ad_alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(inout) :: ad_pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k

      real(r8) :: SedDen, cff, cff1, cff2
      real(r8) :: ad_SedDen, ad_cff, ad_cff1
      real(r8) :: adfac, adfac1

#  ifdef VAR_RHO_2D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhoA1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhoS1
#  endif

#  include "set_bounds.h"
!
!=======================================================================
!  Adjoint linear equation of state.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
#   ifdef BV_FREQUENCY_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_bvf)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_bvf)
#   endif
#   ifdef VAR_RHO_2D
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_rhoA, tl_rhoS)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_rhoA, ad_rhoS)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES_NOT_YET
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_alpha, tl_beta)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_alpha, ad_beta)
#    ifdef LMD_DDMIX_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_alfaobeta)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_alfaobeta)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_pden)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_pden)
#   endif
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_rho)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_rho)
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#  ifdef BV_FREQUENCY_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          tl_bvf)
!>
        CALL ad_exchange_w3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                             ad_bvf)
#  endif
#  ifdef VAR_RHO_2D
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_rhoS)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_rhoS)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_rhoA)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_rhoA)
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_beta)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_beta)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_alpha)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_alpha)
#   ifdef LMD_DDMIX_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          tl_alfaobeta)
!>
        CALL ad_exchange_w3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                             ad_alfaobeta)
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_pden)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_pden)
#  endif
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_rho)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_rho)
      END IF
!
      DO j=JstrT,JendT

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute adjoint thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------

#   ifdef LMD_DDMIX_NOT_YET
!
!  Compute ratio of thermal expansion and saline contraction
!  coefficients.
!
        IF (Scoef(ng).eq.0.0_r8) THEN
          cff=1.0_r8
        ELSE
          cff=1.0_r8/Scoef(ng)
        END IF
        DO k=1,N(ng)
          DO i=IstrT,IendT
!>          tl_alfaobeta(i,j,k)=0.0_r8
!>
            ad_alfaobeta(i,j,k)=0.0_r8
          END DO
        END DO
#   endif
!
!  Set thermal expansion and saline contraction coefficients.
!
        DO i=IstrT,IendT
#   ifdef SALINITY
!>        tl_beta(i,j)=0.0_r8
!>
          ad_beta(i,j)=0.0_r8
#   else
!>        beta(i,j)=0.0_r8
!>
          ad_beta(i,j)=0.0_r8
#   endif
!>        tl_alpha(i,j)=0.0_r8
!>
          ad_alpha(i,j)=0.0_r8
        END DO
#  endif

#  ifdef BV_FREQUENCY_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points.
!-----------------------------------------------------------------------
!
        ad_cff=0.0_r8
        DO k=1,N(ng)
          DO i=IstrT,IendT
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
!>          tl_bvf(i,j,k)=-gorho0*                                      &
!>   &                    (cff*(tl_rho(i,j,k+1)-tl_rho(i,j,k))+         &
!>   &                     tl_cff*(rho(i,j,k+1)-rho(i,j,k)))
!>
            adfac=-gorho0*ad_bvf(i,j,k)
            adfac1=adfac*cff
            ad_rho(i,j,k+1)=ad_rho(i,j,k+1)+adfac1
            ad_rho(i,j,k  )=ad_rho(i,j,k  )-adfac1
            ad_cff=ad_cff+(rho(i,j,k+1)-rho(i,j,k))*adfac
            ad_bvf(i,j,k)=0.0_r8
!>          tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
!>
            adfac=-cff*cff*ad_cff
            ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)+adfac
            ad_z_r(i,j,k  )=ad_z_r(i,j,k  )-adfac
            ad_cff=0.0_r8
          END DO
        END DO
#  endif

#  ifdef VAR_RHO_2D
!
!---------------------------------------------------------------------
!  Compute adjoint vertical averaged density (ad_rhoA) and adjoint
!  density perturbation (ad_rhoS) used in adjoint barotropic pressure
!  gradient.
!---------------------------------------------------------------------
!
!  Compute intermediate BASIC STATE rhoS1 and rhoA1.
!
        DO i=IstrT,IendT
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoS1(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA1(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=rho(i,j,k)*Hz(i,j,k)
            rhoS1(i,j)=rhoS1(i,j)+Hz(i,j,k)*(rhoA1(i,j)+0.5_r8*cff1)
            rhoA1(i,j)=rhoA1(i,j)+cff1
          END DO
        END DO
!
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
!>        tl_rhoS(i,j)=2.0_r8*cff2*                                     &
!>   &                 cff1*(2.0_r8*tl_cff1*rhoS1(i,j)+                 &
!>   &                       cff1*tl_rhoS(i,j))
!>
          adfac=2.0_r8*cff2*cff1*ad_rhoS(i,j)
          ad_cff1=2.0_r8*rhoS1(i,j)*adfac
          ad_rhoS(i,j)=cff1*adfac
!>        tl_rhoA(i,j)=cff2*(tl_cff1*rhoA1(i,j)+cff1*tl_rhoA(i,j))
!>
          adfac=cff2*ad_rhoA(i,j)
          ad_cff1=ad_cff1+rhoA1(i,j)*adfac
          ad_rhoA(i,j)=cff1*adfac
!>        tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
!>
          adfac=-cff1*cff1*ad_cff1
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac
          ad_z_w(i,j,0    )=ad_z_w(i,j,0    )-adfac
          ad_cff1=0.0_r8
        END DO
!
!  Compute appropriate intermediate BASIC STATE "rhoA1".
!
        DO i=IstrT,IendT
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoA1(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=rho(i,j,k)*Hz(i,j,k)
!>          tl_rhoA(i,j)=tl_rhoA(i,j)+tl_cff1
!>
            ad_cff1=ad_rhoA(i,j)
!>          tl_rhoS(i,j)=tl_rhoS(i,j)+                                  &
!>   &                   tl_Hz(i,j,k)*(rhoA1(i,j)+0.5_r8*cff1)+         &
!>   &                   Hz(i,j,k)*(tl_rhoA(i,j)+0.5_r8*tl_cff1)
!>
            adfac=Hz(i,j,k)*ad_rhoS(i,j)
            ad_rhoA(i,j)=ad_rhoA(i,j)+adfac
            ad_cff1=ad_cff1+0.5_r8*adfac
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                  &
     &                   (rhoA1(i,j)+0.5_r8*cff1)*ad_rhoS(i,j)
!>          tl_cff1=tl_rho(i,j,k)*Hz(i,j,k)+                            &
!>   &              rho(i,j,k)*tl_Hz(i,j,k)
!>
            ad_rho(i,j,k)=ad_rho(i,j,k)+Hz(i,j,k)*ad_cff1
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+rho(i,j,k)*ad_cff1
            ad_cff1=0.0_r8
            rhoA1(i,j)=rhoA1(i,j)+cff1
          END DO
        END DO
        DO i=IstrT,IendT
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
!>        tl_rhoA(i,j)=tl_cff1
!>
          ad_cff1=ad_rhoA(i,j)
          ad_rhoA(i,j)=0.0_r8
!>        tl_rhoS(i,j)=0.5_r8*(tl_cff1*Hz(i,j,N(ng))+
!>   &                         cff1*tl_Hz(i,j,N(ng)))
!>
          adfac=0.5_r8*ad_rhoS(i,j)
          ad_cff1=ad_cff1+Hz(i,j,N(ng))*adfac
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+cff1*adfac
          ad_rhoS(i,j)=0.0_r8
!>        tl_cff1=tl_rho(i,j,N(ng))*Hz(i,j,N(ng))+
!>   &            rho(i,j,N(ng))*tl_Hz(i,j,N(ng))
!>
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))+Hz(i,j,N(ng))*ad_cff1
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+rho(i,j,N(ng))*ad_cff1
          ad_cff1=0.0_r8
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Compute adjoint "in situ" density anomaly (kg/m3 - 1000) using linear
!  equation of state.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>          tl_pden(i,j,k)=tl_rho(i,j,k)
!>
            ad_rho(i,j,k)=ad_rho(i,j,k)+ad_pden(i,j,k)
            ad_pden(i,j,k)=0.0
#  endif
#  ifdef MASKING
!>          tl_rho(i,j,k)=tl_rho(i,j,k)*rmask(i,j)
!>
            ad_rho(i,j,k)=ad_rho(i,j,k)*rmask(i,j)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_DENS_NOT_YET
!>          tl_rho(i,j,k)=tl_rho(i,j,k)+tl_SedDen
!>
            ad_SedDen=ad_SedDen+tl_rho(i,j,k)
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
!>            tl_SedDen=tl_SedDen+                                      &
!>   &                  (tl_t(i,j,k,nrhs,itrc)*                         &
!>   &                   (Srho(ised,ng)-rho(i,j,k))-                    &
!>   &                   t(i,j,k,nrhs,itrc)*                            &
!>   &                   tl_rho(i,j,k))*cff1
!>
              adfac=cff1*ad_SedDen
              ad_rho(i,j,k)=ad_rho(i,j,k)-                              &
     &                      t(i,j,k,nrhs,itrc)*adfac
              tl_t(i,j,k,nrhs,itrc)=tl_t(i,j,k,nrhs,itrc)+              &
     &                              (Srho(ised,ng)-rho(i,j,k))*adfac
            END DO
!>          tl_SedDen=0.0_r8
!>
            ad_SedDen=0.0_r8
#  endif
#  ifdef SALINITY
!>          tl_rho(i,j,k)=tl_rho(i,j,k)+                                &
!>   &                    R0(ng)*Scoef(ng)*t(i,j,k,nrhs,isalt)
!>
            ad_t(i,j,k,nrhs,isalt)=ad_t(i,j,k,nrhs,isalt)+              &
     &                             R0(ng)*Scoef(ng)*ad_rho(i,j,k)
#  endif
!>          tl_rho(i,j,k)=-R0(ng)*Tcoef*tl_t(i,j,k,nrhs,itemp)
!>
            ad_t(i,j,k,nrhs,itemp)=ad_t(i,j,k,nrhs,itemp)-              &
     &                             R0(ng)*Tcoef(ng)*ad_rho(i,j,k)
            ad_rho(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE ad_rho_eos_tile
# endif
#endif
      END MODULE ad_rho_eos_mod
