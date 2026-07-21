      MODULE biology_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group            VIMS Colleagues   !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.md                                Katja Fennel   !
!=======================================================================
!                                                                      !
!  This routine computes the  biological sources and sinks terms for   !
!  the Estuarine Carbon Biogeochemistry model (Feng et al., 2015, and  !
!  follow-ups).                                                        !
!                                                                      !
!  This model is loosely based on the model by Fasham et al. (1990)    !
!  but it differs in many respects.  The detailed equations of the     !
!  nitrogen cycling component  are given in  Fennel et al. (2006),     !
!  Feng et al. (2015), and multiple follow-ups listed below.           !
!  Nitrogen is the fundamental elemental currency in this model.       !
!                                                                      !
!  Available options:                                                  !
!                                                                      !
!  - "pCO2_RZ":                                                        !
!                                                                      !
!    Computes the equilibrium partial pressure of CO2 from the         !
!    routines by Zeebe and Wolf-Gladrow (2001), while the OCMIP        !
!    standard routines are the default. There are two different        !
!    ways of treating alkalinity. It can be treated diagnostically     !
!    (default), in this case alkalinity acts like a passive tracer     !
!    that is not affected by changes in the concentration of nitrate   !
!    or ammonium. However, if the option "TALK_NONCONSERV" is used,    !
!    the alkalinity will be affected by sources and sinks in nitrate.  !
!                                                                      !
!  - "OCMIP_OXYGEN_SC":                                                !
!                                                                      !
!    The Schmidt number of oxygen in seawater will be computed using   !
!    the formulation  proposed by Keeling et al. (1998). Otherwise,    !
!    the Wanninkhof formula will be used.                              !
!                                                                      !
!  - "PO4":                                                            !
!                                                                      !
!    The phytoplankton growth can be limited by either nitrogen or     !
!    phosphorus. This option was introduced in Laurent et al. (2012).  !
!                                                                      !
!  - "RW14_OXYGEN_SC" and/or  "RW14_CO2_SC":                           !
!                                                                      !
!    The model use a Wanninkhof (2014) air-sea flux parameterization.  !
!                                                                      !
!  - "PCO2AIR_DATA":                                                   !
!                                                                      !
!    The atmospheric pCO2 uses Eq.1 from St-Laurent et al. (2020).     !
!                                                                      !
!  - "PCO2AIR_SECULAR":                                                !
!                                                                      !
!    An alternative time-dependent atmospheric pCO2 evolution. If      !
!    neither options for pCO2 is activated, the atmospheric pCO2       !
!    is constant.                                                      !
!                                                                      !
!  The "ECB" model was forked from "BIO_FENNEL" some time after        !
!  Fennel et al. (2008). Multiple processes specific to estuaries      !
!  were added in Druon et al. (2010), Feng et al. (2015), and in       !
!  multiple follow-ups listed below.                                   !
!                                                                      !
!  The "ECB" code attempts to preserve, as muchs as possible, the      !
!  structure and content of github.com/myroms "BIO_FENNEL" code, in    !
!  order to facilitate updates with the ROMS repository. Thus, "ECB"   !
!  and "BIO_FENNEL" share a large portion of their code, as it become  !
!  apparent when  comparing them side-by-side. Nevertheless, "ECB"     !
!  and "BIO_FENNEL" remain distinct code bases and should be           !
!  considered as such.                                                 !
!                                                                      !
!  The options "CARBON","OXYGEN", "BIO_SEDIMENT", "DENITRIFICATION"    !
!  should always be activated when using "ECB". The lines of code      !
!  below implicitly expects these options to be activated in the       !
!  header file, which allowed for a considerable reduction in the      !
!  number of CPP directives compared to "BIO_FENNEL".                  !
!                                                                      !
!**********************************************************************!
!                                                                      !
!  References relevant for "BIO_FENNEL" and "ECB":                     !
!                                                                      !
!    Fennel, K., Wilkin, J., Levin, J., Moisan, J., O^Reilly, J.,      !
!      Haidvogel, D., 2006: Nitrogen cycling in the Mid Atlantic       !
!      Bight and implications for the North Atlantic nitrogen          !
!      budget: Results from a three-dimensional model, Global          !
!      Biogeochemical Cycles, 20, GB3007, doi:10.1029/2005GB002456.    !
!                                                                      !
!    Fennel, K., Wilkin, J., Previdi, M., Najjar, R. 2008:             !
!      Denitrification effects on air-sea CO2 flux in the coastal      !
!      ocean: Simulations for the Northwest North Atlantic,            !
!      Geophys. Res. Letters 35, L24608, doi:10.1029/2008GL036147.     !
!                                                                      !
!    Fennel, K., Hu, J., Laurent, A., Marta-Almeida, M., Hetland, R.   !
!      2013: Sensitivity of Hypoxia Predictions for the Northern Gulf  !
!      of Mexico to Sediment Oxygen Consumption and Model Nesting, J.  !
!      Geophys. Res. Ocean 118 (2), 990-1002, doi:10.1002/jgrc.20077.  !
!                                                                      !
!    Keeling, R.F., Stephens, B.B., Najjar, R.G., Doney, S.C.,         !
!      Archer, D., and Heimann, M., 1998: Seasonal variations in the   !
!      atmospheric O2/N2 ratio in relation to the kinetics of air-sea  !
!      gas exchange, Global Biogeochem. Cycles, 12, 141-163.           !
!                                                                      !
!    Laurent, A., Fennel, K., Hu, J., Hetland, R. 2012: Simulating     !
!      the Effects of Phosphorus Limitation in the Mississippi and     !
!      Atchafalaya River Plumes, Biogeosciences, 9 (11), 4707-4723,    !
!      doi:10.5194/bg-9-4707-2012.                                     !
!                                                                      !
!    Yu, L., Fennel, K., Laurent, A., Murrell, M. C., Lehrter, J. C.   !
!      2015: Numerical Analysis of the Primary Processes Controlling   !
!      Oxygen Dynamics on the Louisiana Shelf, Biogeosciences, 12 (7), !
!      2063-2076, doi:10.5194/bg-12-2063-2015.                         !
!                                                                      !
!    Wanninkhof, R., 2014: Relationship between Wind Speed and Gas     !
!      Exchange over the Ocean Revisited, Limnol. Oceanogr. Methods    !
!      12 (6), 351-362, doi:10.4319/lom.2014.12.351.                   !
!                                                                      !
!**********************************************************************!
!                                                                      !
!  References specific to "ECB" (Estuarine Carbon Biogeochemistry):    !
!                                                                      !
!    Druon, J.N., Mannino, A., Signorini, S., McClain, Friedrich, M.,  !
!      Wilkin, J., and Fennel, K., 2010: Modeling the dynamics and     !
!      export of dissolved organic matter in the Northeastern U.S.     !
!      continental shelf, Estuarine, Coastal, and Shelf Science, 88,   !
!      488-507, doi:10.1016/j.ecss.2010.05.010.                        !
!                                                                      !
!    Feng,Y., Fiedricks, M., Wilkin, J., Tian, H., Yang, Q., Hofmann,  !
!      A.E., Wiggert, J.D., Hood, R.R., 2015: Chesapeake Bay nitrogen  !
!      fluxes derived from a land-estuarine ocean biogeochemical       !
!      modeling system: Model description, evaluation, and nitrogen    !
!      budgets, J. Geophys. Res. Biogeosci., 120, 1666-1695,           !
!      doi:10.1002/2015JG002931.                                       !
!                                                                      !
!    Da, F., Friedrichs, M.A.M., St-Laurent, P. (2018). Impacts of     !
!      atmospheric nitrogen deposition and coastal nitrogen fluxes on  !
!      oxygen concentrations in Chesapeake Bay, J. Geophys. Res.       !
!      Oceans, 123, 5004-5025. doi:10.1029/2018JC014009.               !
!                                                                      !
!    St-Laurent, P., Friedrichs, M.A.M., Najjar, R.G., Shadwick, E.H., !
!      Tian, H., and Yao, Y., 2020: Relatice impact of global changes  !
!      on the inorganic carbon balance of the Chesapeake Bay,          !
!      Biogeosciences, 17, 3779-3796, doi:10.5194/bg-17-3779-2020.     !
!                                                                      !
!    Turner, J.S., St-Laurent, P., Friedrichs, M.A.M., Friedrichs,     !
!      C.T., 2021: Effects of reduced shoreline erosion on Chesapeake  !
!      Bay water clarity, J. Sci. T. Env., 769, 145157,                !
!      doi:10.1016/j.scitotenv.2021.145157.                            !
!                                                                      !
!    Da, F., Friedrichs, M.A.M., St-Laurent, P., Shadwick, E.H.,       !
!      Najjar, R.G., Hinson, K.E., 2021: Mechanisms driving decadal    !
!      changes in the carbonate system of a coastal plain estuary,     !
!      J. Geophys. Res. Oceans, 126, e2021JC017239                     !
!      doi:10.1029/2021jc017239.                                       !
!                                                                      !
!    Da, F., Friedrichs, M.A.M., St-Laurent, P., Najjar, R.G.,         !
!      Shadwick, E.H., Stets, E.G., 2024: Influence of Rivers, Tides,  !
!      and Tidal Wetlands on Estuarine Carbonate System Dynamics,      !
!      Estuaries Coast., 47, 2283-2305, doi:10.1007/s12237-024-01421-z !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_biology
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE dateclock_mod, ONLY : caldate
#if defined SEDIMENT && defined SANDS_PROXY_ISS
      USE mod_sediment,  ONLY : isand
#endif
!
      implicit none
!
      PUBLIC  :: biology
      PRIVATE :: ecb_tile
#ifdef pCO2_RZ
      PRIVATE :: pCO2_water_RZ
#else
      PRIVATE :: pCO2_water
#endif
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE biology (ng, tile)
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
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=MyFile
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15, __LINE__, MyFile)
#endif
      CALL ecb_tile (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj, N(ng), NT(ng),                 &
     &               IminS, ImaxS, JminS, JmaxS,                        &
     &               nstp(ng), nnew(ng),                                &
#ifdef MASKING
     &               GRID(ng) % rmask,                                  &
# ifdef WET_DRY
     &               GRID(ng) % rmask_wet,                              &
#  ifdef DIAGNOSTICS_BIO
     &               GRID(ng) % rmask_full,                             &
#  endif
# endif
#endif
     &               GRID(ng) % Hz,                                     &
     &               GRID(ng) % z_r,                                    &
     &               GRID(ng) % z_w,                                    &
     &               FORCES(ng) % srflx,                                &
#ifdef BULK_FLUXES
     &               FORCES(ng) % Uwind,                                &
     &               FORCES(ng) % Vwind,                                &
#else
     &               FORCES(ng) % sustr,                                &
     &               FORCES(ng) % svstr,                                &
#endif
     &               FORCES(ng) % bustr,                                &
     &               FORCES(ng) % bvstr,                                &
     &               OCEAN(ng) % pH,                                    &
#ifdef DIAGNOSTICS_BIO
     &               DIAGS(ng) % DiaBio2d,                              &
     &               DIAGS(ng) % DiaBio3d,                              &
#endif
     &               OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE biology
!
!***********************************************************************
      SUBROUTINE ecb_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, UBk, UBt,                &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
#ifdef MASKING
     &                     rmask,                                       &
# if defined WET_DRY
     &                     rmask_wet,                                   &
#  ifdef DIAGNOSTICS_BIO
     &                     rmask_full,                                  &
#  endif
# endif
#endif
     &                     Hz, z_r, z_w, srflx,                         &
#ifdef BULK_FLUXES
     &                     Uwind, Vwind,                                &
#else
     &                     sustr, svstr,                                &
#endif
                           bustr, bvstr,                                &
     &                     pH,                                          &
#ifdef DIAGNOSTICS_BIO
     &                     DiaBio2d, DiaBio3d,                          &
#endif
     &                     t)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#   ifdef DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
#   endif
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
# else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
# endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#   ifdef DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
#   endif
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 6

      integer :: Iter, i, ibio, ised, isink, itrc, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

      real(r8) :: u10squ

#if defined OCMIP_OXYGEN_SC
!
! Alternative formulation for Schmidt number coefficients (Sc will be
! slightly smaller up to about 35C) using the formulation proposed by
! Keeling et al. (1998, Global Biogeochem. Cycles, 12, 141-163).
!
      real(r8), parameter :: A_O2 = 1638.0_r8
      real(r8), parameter :: B_O2 = 81.83_r8
      real(r8), parameter :: C_O2 = 1.483_r8
      real(r8), parameter :: D_O2 = 0.008004_r8
      real(r8), parameter :: E_O2 = 0.0_r8

#elif defined RW14_OXYGEN_SC
!
! Alternative formulation for Schmidt number coefficients using the
! formulation of Wanninkhof (2014, L and O Methods, 12,351-362).
!
      real(r8), parameter :: A_O2 = 1920.4_r8
      real(r8), parameter :: B_O2 = 135.6_r8
      real(r8), parameter :: C_O2 = 5.2122_r8
      real(r8), parameter :: D_O2 = 0.10939_r8
      real(r8), parameter :: E_O2 = 0.00093777_r8

#else
!
! Schmidt number coefficients using the formulation of
! Wanninkhof (1992).
!
      real(r8), parameter :: A_O2 = 1953.4_r8
      real(r8), parameter :: B_O2 = 128.0_r8
      real(r8), parameter :: C_O2 = 3.9918_r8
      real(r8), parameter :: D_O2 = 0.050091_r8
      real(r8), parameter :: E_O2 = 0.0_r8
#endif
      real(r8), parameter :: OA0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: OA1 = 3.22014_r8       ! saturation
      real(r8), parameter :: OA2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: OA3 = 4.94457_r8
      real(r8), parameter :: OA4 =-0.256847_r8
      real(r8), parameter :: OA5 = 3.88767_r8
      real(r8), parameter :: OB0 =-0.00624523_r8
      real(r8), parameter :: OB1 =-0.00737614_r8
      real(r8), parameter :: OB2 =-0.0103410_r8
      real(r8), parameter :: OB3 =-0.00817083_r8
      real(r8), parameter :: OC0 =-0.000000488682_r8
      real(r8), parameter :: rOxNO3= 8.625_r8       ! 138/16
      real(r8), parameter :: rOxNH4= 6.625_r8       ! 106/16

! Parameters K_NTR and K_DNF are used to partition processes between
! 'aerobic' and 'anaerobic' through the functions f_NTR and f_DNF
! (Feng et al., 2015). This partitioning only makes sense if
! 0<=f_NTR<=1, 0<=f_DNF<=1, and f_NTR+f_DNF<=1. An easy way to
! enforce this constraint is to have K_DNF=K_NTR (as was the case
! in Feng et al., 2015; Table-A4).
!
      real(r8), parameter :: K_NTR = 1.0_r8      ! Feng et al., 2015
      real(r8), parameter :: K_DNF = K_NTR       ! Enforce the equality
      real(r8), parameter :: K_WNO3 = 0.5_r8
      real(r8), parameter :: K_BO2 = 26.5_r8
      real(r8) :: l2mol = 1000.0_r8/22.3916_r8   ! liter to mol
      integer :: year
      integer, parameter :: DoNewton = 0         ! pCO2 solver
!
#if defined RW14_CO2_SC
      real(r8), parameter :: A_CO2 = 2116.8_r8      ! Schmidt number
      real(r8), parameter :: B_CO2 = 136.25_r8      ! transfer coeff
      real(r8), parameter :: C_CO2 = 4.7353_r8      ! according to
      real(r8), parameter :: D_CO2 = 0.092307_r8    ! Wanninkhof (2014)
      real(r8), parameter :: E_CO2 = 0.0007555_r8
#else
      real(r8), parameter :: A_CO2 = 2073.1_r8      ! Schmidt
      real(r8), parameter :: B_CO2 = 125.62_r8      ! number
      real(r8), parameter :: C_CO2 = 3.6276_r8      ! transfer
      real(r8), parameter :: D_CO2 = 0.043219_r8    ! coefficients
      real(r8), parameter :: E_CO2 = 0.0_r8
#endif
!
      real(r8), parameter :: A1 = -60.2409_r8       ! surface
      real(r8), parameter :: A2 = 93.4517_r8        ! CO2
      real(r8), parameter :: A3 = 23.3585_r8        ! solubility
      real(r8), parameter :: B1 = 0.023517_r8       ! coefficients
      real(r8), parameter :: B2 = -0.023656_r8
      real(r8), parameter :: B3 = 0.0047036_r8
!
      real(r8) :: pmonth                         ! months since Jan 1951
      real(r8) :: pCO2air_secular
      real(dp) :: yday

      real(r8), parameter :: pi2 = 6.2831853071796_r8

      real(r8), parameter :: D0 = 282.6_r8          ! coefficients
      real(r8), parameter :: D1 = 0.125_r8          ! to calculate
      real(r8), parameter :: D2 =-7.18_r8           ! secular trend in
      real(r8), parameter :: D3 = 0.86_r8           ! atmospheric pCO2
      real(r8), parameter :: D4 =-0.99_r8
      real(r8), parameter :: D5 = 0.28_r8
      real(r8), parameter :: D6 =-0.80_r8
      real(r8), parameter :: D7 = 0.06_r8

      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
#ifdef PO4
      real(r8), parameter :: MinVal = 1.0e-6_r8

      real(r8) :: L_PO4, LMIN, mu
#endif
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: cff7, cff8
      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8) :: total_N

#ifdef DIAGNOSTICS_BIO
      real(r8) :: fiter
#endif
!
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA

      real(r8) :: C_Flux_RemineL, C_Flux_RemineS, C_Flux_RemineR
      real(r8) :: CO2_Flux, CO2_sol, SchmidtN, TempK

      real(r8) :: f_NTR, f_DNF
      real(r8) :: LBO2, O2satB
      real(r8) :: N_Flux_Assim
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_Egest
      real(r8) :: N_Flux_NewProd, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_RemineL, N_Flux_RemineS, N_Flux_RemineR
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo
      real(r8) :: TSS, ReSuspR, Cbe, Ustarb
      real(r8) :: C_exc_Up, qDOC, rDOC, CNbur
      real(r8) :: g_gmax, ZooAE_N_g, qDON, rDON
      real(r8) :: N_Flux_Sloppy_lDON, N_Flux_Sloppy_slDON
      real(r8), dimension(Nsink) :: Wbio
!
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: pCO2

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

#include "set_bounds.h"

#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)

#ifdef DIAGNOSTICS_BIO
!
!  A factor to account for the number of iterations in accumulating
!  diagnostic rate variables.
!
      fiter=1.0_r8/REAL(BioIter(ng),r8)
#endif
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt
      idsink(2)=iChlo
      idsink(3)=iSDeN
      idsink(4)=iLDeN
      idsink(5)=iSDeC
      idsink(6)=iLDeC
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! phytoplankton
      Wbio(2)=wPhy(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus
      Wbio(5)=wSDet(ng)               ! small Carbon-detritus
      Wbio(6)=wLDet(ng)               ! large Carbon-detritus
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8, t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO

#if defined SEDIMENT && defined SANDS_PROXY_ISS
!
        DO ised=1,NNS
          itrc=isand(ised)     ! suspended noncohesive sediment tracers
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,itrc)=MAX(t(i,j,k,nstp,itrc), 0.0_r8)
            END DO
          END DO
        END DO
#endif
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio_old(i,k,iTIC_)=MIN(Bio_old(i,k,iTIC_),3000.0_r8)
            Bio_old(i,k,iTIC_)=MAX(Bio_old(i,k,iTIC_),  1.0_r8)
            Bio(i,k,iTIC_)=Bio_old(i,k,iTIC_)
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in the sense that the sum of concentration
!  of all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  substances    within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of substances    above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
#if defined SEDIMENT && defined SANDS_PROXY_ISS
                TSS=(Bio(i,k,iZoop)*ZooCN(ng) +                         &
     &               Bio(i,k,iPhyt)*PhyCN(ng) +                         &
     &               Bio(i,k,iSDeC)           +                         &
     &               Bio(i,k,iLDeC))*12.0_r8*1.e-3_r8     ! POC, g-C/m3
!
!  Turner et al., 2021: The total suspended solids (TSS) is computed
!  as the sum of volatil (VSS) and fixed (FSS) components, TSS=VSS+FSS,
!  where VSS=POC*2.9 and FSS=0.35*VSS+ISS. Here, suspended noncohesive
!  sediments (sand size classes) are used as a proxy for the inorganic
!  suspended solids (ISS).
!
!  Convert POC to VSS (g-C/m3 to g/m3) by using the factor 2.9 of
!  Cerco and Noel (WQSTM, 2017) and then add 0.35*VSS to represent
!  the organic portion of FSS.
!
                TSS=TSS*2.9_r8*1.35_r8       ! 1.35xVSS (g/m3)
!
!  Add ISS contribution. Sand tracers converted from kg/m3 to g/m3.
!
                DO ised=1,NNS
                  itrc=isand(ised)     ! suspended noncohesive sediment
                  TSS=TSS+Bio(i,k,itrc)*1.e3_r8
                END DO
                AttFac=rkdTSS1(ng)*TSS+                                 &
     &                 rkdS1(ng)*Bio(i,k,isalt)
!
!  Enforce a lower bound (salinity term could produce Att < 0).
!
                Att=MAX(rkd1(ng)+AttFac, 0.3_r8)*                       &
     &              (z_w(i,j,k)-z_w(i,j,k-1))
#else
!
!  CDOM contribution: Eq.4 of Feng et al., 2015.
!
                AttFac=MAX(Bio(i,k,iDOC_)+Bio(i,k,irDOC)-70.819_r8,     &
     &                     0.0_r8)*                                     &
     &                 3.786e-3_r8
                Att=(AttSW(ng)+                                         &
     &               AttChl(ng)*Bio(i,k,iChlo)+                         &
     &               AttFac)*                                           &
     &               (z_w(i,j,k)-z_w(i,j,k-1))
#endif
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8
                Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),      &
     &                    Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate with Q10 based on
!  Lomas et al., 2002.
!
                Vp=max(Vp0(ng)*exp(0.08065_r8 * Bio(i,k,itemp)),4.00_r8)
                fac1=PAR*PhyIS(ng)
                Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                t_PPmax=Epp*fac1
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
                cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                inhNH4=1.0_r8/(1.0_r8+cff1)
                L_NH4=cff1/(1.0_r8+cff1)
                L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                LTOT=L_NO3+L_NH4
#ifdef PO4
                cff3=Bio(i,k,iPO4_)*K_PO4(ng)
                L_PO4=cff3/(1.0_r8+cff3)
                LMIN=MIN(LTOT,L_PO4)
#endif
!
!  Nitrate and ammonium uptake by Phytoplankton.
!  Phytoplankton exudation of semilabile DON to DON
!  and phytoplankton exudation of labile DON to NH4
!
#ifdef PO4
                mu=dtdays*t_PPmax*LMIN
                cff4=mu*Bio(i,k,iPhyt)*L_NO3/                           &
     &               MAX(MinVal,LTOT)/MAX(MinVal,Bio(i,k,iNO3_))
                cff5=mu*Bio(i,k,iPhyt)*L_NH4/                           &
     &               MAX(MinVal,LTOT)/MAX(MinVal,Bio(i,k,iNH4_))
                cff6=R_P2N(ng)*mu*Bio(i,k,iPhyt)/                       &
     &               MAX(MinVal,Bio(i,k,iPO4_))
#else
                fac1=dtdays*t_PPmax
                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*Bio(i,k,iPhyt)
                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
#endif
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
#ifdef PO4
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff6)
#endif
                N_Flux_NewProd=Bio(i,k,iNO3_)*cff4
                N_Flux_RegProd=Bio(i,k,iNH4_)*cff5
!
                f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + K_NTR)
                f_DNF = K_DNF/(Bio(i,k,iOxyg) + K_DNF)                  &
     &                * Bio(i,k,iNO3_) / (Bio(i,k,iNO3_) + K_WNO3)
!
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         (1.0_r8-EsDON(ng)-                       &
     &                             (f_NTR+f_DNF)*ElDON(ng))*            &
     &                        (N_Flux_NewProd+N_Flux_RegProd)
!
                Bio(i,k,iChlo)=Bio(i,k,iChlo)+                          &
#ifdef PO4
     &                         (dtdays*t_PPmax*t_PPmax*LMIN*LMIN*       &
#else
     &                         (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*       &
#endif
     &                          Chl2C_m(ng)*Bio(i,k,iChlo))/            &
     &                         (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iPPro)=DiaBio3d(i,j,k,iPPro)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  &
# endif
     &                                (N_Flux_NewProd+N_Flux_RegProd)*  &
     &                                fiter
                DiaBio3d(i,j,k,iNO3u)=DiaBio3d(i,j,k,iNO3u)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  &
# endif
     &                                N_Flux_NewProd*fiter
#endif
!
!  Phytoplankton exudation to semilabileDON and to `labile DON' (NH4)
!
                Bio(i,k,iDON_)=Bio(i,k,iDON_)+EsDON(ng)*                &
     &                        (N_Flux_NewProd+N_Flux_RegProd)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+(f_NTR+f_DNF)*ElDON(ng)*  &
     &                        (N_Flux_NewProd+N_Flux_RegProd)
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+(f_NTR+f_DNF)*ElDON(ng)*  &
     &                        (N_Flux_NewProd+N_Flux_RegProd)
!
!  Remineralization of labile DOM consumes oxygen
!
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
     &                         rOxNO3*N_Flux_NewProd+                   &
     &                         rOxNH4*N_Flux_RegProd-                   &
     &                         rOxNH4*f_NTR*ElDON(ng)*                  &
     &                         (N_Flux_NewProd+N_Flux_RegProd)
                Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
!
!  Total inorganic carbon (CO2) uptake during phytoplankton growth.
!  (Ensure that phytoplankton uptake cannot exceed TIC inventory).
!
                cff1=PhyCN(ng)*(N_Flux_NewProd+N_Flux_RegProd)
                cff1=MIN(cff1, Bio(i,k,iTIC_)-eps)
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-cff1
!
!  The carbon excess uptake (C_exc_Up) is equal to the phytoplankton
!  exudation of total DOC. The C excess uptake is defined as
!  gammaC * [PP(limited by light)-PP(limited by light+nutrient)].
!  Only a fraction of it is semilabile (slCexc, directed to DOC),
!  the rest is labile and is directed back to TIC.
!  (Ensure that C_exc_Up uptake cannot exceed the TIC inventory).
!
                C_exc_Up=gammaC(ng)*                                    &
     &                   PhyCN(ng)*dtdays*t_PPmax*Bio(i,k,iPhyt)*       &
     &                   (1.0_r8-LTOT)
                C_exc_Up=MIN(C_exc_Up,                                  &
     &                       MAX(Bio(i,k,iTIC_), 0.0_r8)/slCexc(ng))
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-C_exc_Up*slCexc(ng)
                Bio(i,k,iDOC_)=Bio(i,k,iDOC_)+C_exc_Up*slCexc(ng)
!
!  The oxygen produced by the synthesis of carbohydrates has a one to
!  one mole ratio with TIC: CO2 + H2O + energy -> CH2O + O2
!
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+C_exc_Up
!
!  Nutrient-based semilabile DOC exudation calculated as a C:N ratio
!  of semilabile DON exudation.
!
                Bio(i,k,iDOC_)=Bio(i,k,iDOC_)                           &
     &                         +PhyCN(ng)*EsDON(ng)*                    &
     &                         (N_Flux_NewProd+N_Flux_RegProd)
!
!  Nutrient-based exudation of labile DOC (directed to TIC)
!
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)                           &
     &                         +PhyCN(ng)*ElDON(ng)*                    &
     &                          (f_NTR + f_DNF)*                        &
     &                         (N_Flux_NewProd+N_Flux_RegProd)

#if defined DIAGNOSTICS_BIO
!
!  The total carbon excess uptake includes the labile fraction of the
!  DOC produced (slCexc is the semilabile fraction)
!
                DiaBio3d(i,j,k,iPPCe)=DiaBio3d(i,j,k,iPPCe)+C_exc_Up*   &
# ifdef WET_DRY
     &                                rmask_full( i, j ) *              &
# endif
     &                                fiter
#endif
#ifdef TALK_NONCONSERV
!
!  Account for the uptake of NO3,NH4 on total alkalinity.
!
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_NewProd-           &
     &                         N_Flux_RegProd
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types.
!
! Nitrification is inhibited at low oxygen concentrations using a
! Michaelis-Menten term.
!
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)     ! O2 max
                fac3=MAX(fac2/(K_NTR+fac2),0.0_r8)  ! MM-O2 dependence
                fac1=dtdays*NitriR(ng)*fac3
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iNifx)=DiaBio3d(i,j,k,iNifx)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  &
# endif
     &                                N_Flux_Nitrifi*fiter
#endif
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#if defined TALK_NONCONSERV
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-2.0_r8*N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE
              DO k=N(ng),1,-1
                f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + K_NTR)
                cff3=dtdays*NitriR(ng) * f_NTR
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iNifx)=DiaBio3d(i,j,k,iNifx)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  &
# endif
     &                                N_Flux_Nitrifi*fiter
#endif
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                Bio(i,k,iOxyg)=MAX(0.0_r8,Bio(i,k,iOxyg))
#if defined TALK_NONCONSERV
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-2.0_r8*N_Flux_Nitrifi
#endif
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L and O 38:468-472]
!-----------------------------------------------------------------------
!
          cff2=dtdays*PhyMR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
! Temp-dependent grazing; 0.0875 <=> Q_10 of 2.4 (Lomas et al., 2002).
!
              fac1=dtdays*ZooGR(ng)*EXP(0.0875_r8*Bio(i,k,itemp))
              cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff3=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
              Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to detritus.
!
              N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
              N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+                            &
     &                       N_Flux_Assim
!
! Sloppy feeding related parameterizations are described in Druon et
! al. (2009; NASA/TM-2009214177)
! Ratio of grazing rate to maximum g/gmax - Druon pg. 11
!
              g_gmax=Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                     &
     &               (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
!
! Grazing-based assimilation efficiency (beta) - Druon pg. 12
!
              ZooAE_N_g=(0.75_r8+ZooBM(ng)-0.65_r8*g_gmax)/             &
     &                  (1.0_r8-ZooER(ng)*g_gmax)
!
! Fraction of DON to DON+PON within the phytoplankton cell.
!
              qDON=0.71_r8*g_gmax
!
! A fraction of unassimilated flux goes to large detritus as fecal
! pellets.
!
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+N_Flux_Egest*(1.0_r8-qDON)
!
! Remaining fraction split between DON (semilabile DON) and NH4
! (labile DON).
!
              N_Flux_Sloppy_lDON=N_Flux_Egest*qDON*(1.0_r8-deltN(ng))
              N_Flux_Sloppy_slDON=N_Flux_Egest*qDON*deltN(ng)
              Bio(i,k,iDON_)=Bio(i,k,iDON_)+N_Flux_Sloppy_slDON
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Sloppy_lDON
              Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_Sloppy_lDON
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
              Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_Pmortal
!
! Fraction of DOC to DOC+POC within the phytoplankton cell.
!
              qDOC=qDON
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)+PhyCN(ng)*N_Flux_Egest*     &
     &                      (1.0_r8-qDOC)
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+PhyCN(ng)*N_Flux_Egest*     &
     &                       qDOC*(1.0_r8-deltC(ng))
              Bio(i,k,iDOC_)=Bio(i,k,iDOC_)+PhyCN(ng)*N_Flux_Egest*     &
     &                       qDOC*deltC(ng)
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       PhyCN(ng)*N_Flux_Pmortal
!
! Flux of carbon excess respiration by zooplankton (to TIC)
! due to C:N ratio difference between phytoplankton and zooplankton
!
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       (PhyCN(ng)-ZooCN(ng))*N_Flux_Assim
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-rOxNH4*N_Flux_Sloppy_lDON
              Bio(i,k,iOxyg)=MAX(0.0_r8, Bio(i,k,iOxyg))
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
          cff1=dtdays*ZooBM(ng)
          fac2=dtdays*ZooMR(ng)
          fac3=dtdays*ZooER(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff2=fac2*Bio(i,k,iZoop)
              cff3=fac1*ZooAE_N(ng)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/                            &
     &                       (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
              N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
#ifdef PO4
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+R_P2N(ng)*N_Flux_Zexcret
#endif
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
              N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef PO4
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+R_P2N(ng)*N_Flux_Zmetabo
#endif
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
              Bio(i,k,iOxyg)=MAX( 0._r8, Bio(i,k,iOxyg) )
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       ZooCN(ng)*N_Flux_Zmortal
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       ZooCN(ng)*(N_Flux_Zmetabo+N_Flux_Zexcret)
#ifdef TALK_NONCONSERV
              Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_Zmetabo+             &
     &                       N_Flux_Zexcret
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
              Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)*cff2
              N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
              N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
     &                       N_Flux_CoagP+N_Flux_CoagD
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)+                            &
     &                       PhyCN(ng)*N_Flux_CoagP+                    &
     &                       Bio(i,k,iSDeC)*cff1
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!  DeN solubilization and DON remineralization to NH4.
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              f_NTR=Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+K_NTR)
              f_DNF=K_DNF/(Bio(i,k,iOxyg)+K_DNF)*                       &
     &              Bio(i,k,iNO3_)/(Bio(i,k,iNO3_)+K_WNO3)
!
              cff=dtdays*EXP(0.0875_r8*Bio(i,k,itemp))
!
              cff1=cff*SDeRRN(ng)
              cff2=cff1*deltN(ng)
              cff3=cff1*(1.0_r8-deltN(ng))
!
              cff4=cff*LDeRRN(ng)
              cff5=cff4*deltN(ng)
              cff6=cff4*(1.0_r8-deltN(ng))
!
              rDON=cff*RDeRRN(ng)
!
              fac=(f_NTR+f_DNF)
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)/(1.0_r8+cff2+cff3*fac)
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)/(1.0_r8+cff5+cff6*fac)
              Bio(i,k,iDON_)=Bio(i,k,iDON_)/(1.0_r8+rDON*fac)
!
              N_Flux_RemineS=Bio(i,k,iSDeN)*cff3*fac
              N_Flux_RemineL=Bio(i,k,iLDeN)*cff6*fac
              N_Flux_RemineR=Bio(i,k,iDON_)*rDON*fac
!
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+                            &
     &                       N_Flux_RemineS+N_Flux_RemineL+             &
     &                       N_Flux_RemineR
#ifdef PO4
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                            &
     &                       R_P2N(ng)*(N_Flux_RemineS+                 &
     &                                  N_Flux_RemineL+                 &
     &                                  N_Flux_RemineR)
#endif
              cff7=(Bio(i,k,iDON_)*rDON+                                &
     &              Bio(i,k,iSDeN)*cff3+                                &
     &              Bio(i,k,iLDeN)*cff6)*rOxNH4/                        &
     &             (Bio(i,k,iOxyg)+K_NTR)
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)/(1.0_r8 + cff7)
#if defined TALK_NONCONSERV
              Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+                            &
     &                       N_Flux_RemineS+                            &
     &                       N_Flux_RemineL+                            &
     &                       N_Flux_RemineR
#endif
!
! NO3 lost due to water column denitrification.
! (eta_DNF=84.8/16=5.3, Table A4, Feng et al., 2015).
!
              cff8=(Bio(i,k,iSDeN)*cff3+                                &
     &              Bio(i,k,iLDeN)*cff6+                                &
     &              Bio(i,k,iDON_)*rDON)*5.3_r8/                        &
     &             (Bio(i,k,iOxyg)+ K_DNF )*K_DNF/                      &
     &             (Bio(i,k,iNO3_)+ K_WNO3)
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff8)
              Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+Bio(i,k,iNO3_)*cff8
#ifdef DIAGNOSTICS_BIO
              DiaBio3d(i,j,k,iwdeN)=DiaBio3d(i,j,k,iwdeN)+              &
     &                              Bio(i,k,iNO3_)*cff8*                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              fiter
#endif
!
! N.B. PON solubilization must come last, so that it does not
! inflate the O2 sink associated with DON remineralization nor the
! NO3 sink associated with water-column denitrification.
!
              Bio(i,k,iDON_)=Bio(i,k,iDON_)+                            &
     &                       Bio(i,k,iSDeN)*cff2+                       &
     &                       Bio(i,k,iLDeN)*cff5
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
          cff1=rho0*550.0_r8
#if defined RW14_OXYGEN_SC
          cff2=dtdays*0.251_r8*24.0_r8/100.0_r8
#else
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
#endif
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
#ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
#else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
#endif
            SchmidtN_Ox=A_O2-                                           &
     &                  Bio(i,k,itemp)*(B_O2-                           &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (C_O2-Bio(i,k,itemp)*           &
     &                                  (D_O2-Bio(i,k,itemp)*E_O2)))
            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L and O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            AA=OA0+                                                     &
     &         TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+              &
     &         Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+           &
     &         OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
            O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
            O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
     &                     O2_Flux*Hz_inv(i,k)
#ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iO2fx)=DiaBio2d(i,j,iO2fx)+                    &
# ifdef WET_DRY
     &                          rmask_full(i,j)*                        &
# endif
     &                          O2_Flux*fiter
#endif

          END DO
!
!-----------------------------------------------------------------------
!  Allow different remineralization rates for detrital C and detrital N.
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              f_NTR=Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+K_NTR)
              f_DNF=K_DNF/(Bio(i,k,iOxyg) + K_DNF)*                     &
     &              Bio(i,k,iNO3_)/(Bio(i,k,iNO3_)+K_WNO3)
!
              cff=dtdays*EXP(0.0875_r8*Bio(i,k,itemp))
              cff1=cff*SDeRRC(ng)
              cff2=cff1*deltC(ng)
              cff3=cff1*(1.0_r8-deltC(ng))
!
              cff4=cff*LDeRRC(ng)
              cff5=cff4*deltC(ng)
              cff6=cff4*(1.0_r8-deltC(ng))
!
              rDOC=cff*RDeRRC(ng)
!
              fac=(f_NTR+f_DNF)
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)/(1.0_r8+cff2+cff3*fac)
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)/(1.0_r8+cff5+cff6*fac)
              Bio(i,k,iDOC_)=Bio(i,k,iDOC_)/(1.0_r8+rDOC*fac)
              C_Flux_RemineS=Bio(i,k,iSDeC)*cff3*fac
              C_Flux_RemineL=Bio(i,k,iLDeC)*cff6*fac
              C_Flux_RemineR=Bio(i,k,iDOC_)*rDOC*fac
!
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       C_Flux_RemineS+C_Flux_RemineL+             &
     &                       C_Flux_RemineR
!
!  Fraction deltC of POC is solubilized into semilabile DOC:
!
              Bio(i,k,iDOC_)=Bio(i,k,iDOC_)+                            &
     &                       cff2*Bio(i,k,iSDeC)+                       &
     &                       cff5*Bio(i,k,iLDeC)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Surface CO2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface.
!
          k=N(ng)
#ifdef pCO2_RZ
          CALL pCO2_water_RZ (Istr, Iend, LBi, UBi, LBj, UBj,           &
     &                        IminS, ImaxS, j, DoNewton,                &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        Bio(IminS:,k,itemp), Bio(IminS:,k,isalt), &
     &                        Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk), &
     &                        pH, pCO2)
#else
          CALL pCO2_water (Istr, Iend, LBi, UBi, LBj, UBj,              &
     &                     IminS, ImaxS, j, DoNewton,                   &
# ifdef MASKING
     &                     rmask,                                       &
# endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
     &                     0.0_r8, 0.0_r8, pH, pCO2)
#endif
!
!  Compute surface CO2 gas exchange.
!
          cff1=rho0*550.0_r8
#if defined RW14_CO2_SC
          cff2=dtdays*0.251_r8*24.0_r8/100.0_r8
#else
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
#endif
          DO i=Istr,Iend
!
!  Compute CO2 transfer velocity : u10squared (u10 in m/s)
!
#ifdef BULK_FLUXES
            u10squ=Uwind(i,j)**2+Vwind(i,j)**2
#else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
#endif
            SchmidtN=A_CO2-Bio(i,k,itemp)*(B_CO2-Bio(i,k,itemp)*(C_CO2- &
     &                                           Bio(i,k,itemp)*(D_CO2- &
     &                                           Bio(i,k,itemp)*E_CO2)))
            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN)
!
!  Calculate CO2 solubility [mol/(kg.atm)] using Weiss (1974) formula.
!
            TempK=0.01_r8*(Bio(i,k,itemp)+273.15_r8)
            CO2_sol=EXP(A1+                                             &
     &                  A2/TempK+                                       &
     &                  A3*LOG(TempK)+                                  &
     &                  Bio(i,k,isalt)*(B1+TempK*(B2+B3*TempK)))
!
!  Add in CO2 gas exchange.
!
            CALL caldate (tdays(ng), yy_i=year, yd_dp=yday)
#if defined PCO2AIR_DATA
!
!   Use Eq. 1 from St-Laurent et al. (2020) and redefine day of the year
!   as years since Jan 1, 2001. HGA: Using a year length of 365.2425 days
!   for the Gregorian calendar.
!
            yday=REAL(year,dp)-2001.0_dp+(yday/365.25_dp)
            pCO2air_secular=1.8607_r8*REAL(yday,r8)+                    &
     &                      371.19_r8+0.012447_r8*REAL(yday,r8)**2
            CO2_Flux=cff3*CO2_sol*(pCO2air_secular-pCO2(i))
#elif defined PCO2AIR_SECULAR
            pmonth=year-1951.0_r8+yday/365.0_r8
            pCO2air_secular=D0+D1*pmonth*12.0_r8+                       &
     &                         D2*SIN(pi2*pmonth+D3)+                   &
     &                         D4*SIN(pi2*pmonth+D5)+                   &
     &                         D6*SIN(pi2*pmonth+D7)
            CO2_Flux=cff3*CO2_sol*(pCO2air_secular-pCO2(i))
#else
            CO2_Flux=cff3*CO2_sol*(pCO2air(ng)-pCO2(i))
#endif
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                              &
     &                     CO2_Flux*Hz_inv(i,k)
#ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                    &
# ifdef WET_DRY
     &                          rmask_full(i,j)*                        &
# endif
     &                          CO2_Flux*fiter
            DiaBio2d(i,j,ipCO2)=pCO2(i)
# ifdef WET_DRY
            DiaBio2d(i,j,ipCO2)=DiaBio2d(i,j,ipCO2)*rmask_full(i,j)
# endif
#endif
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=max( Bio(i,k,ibio), 0._r8 )
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
            cff2=4.0_r8/16.0_r8
            cff3=115.0_r8/16.0_r8
            cff4=106.0_r8/16.0_r8
            DO i=Istr,Iend
!
              IF ((ibio.eq.iPhyt).or.                                   &
     &            (ibio.eq.iSDeC).or.                                   &
     &            (ibio.eq.iLDeC).or.                                   &
     &            (ibio.eq.iSDeN).or.                                   &
     &            (ibio.eq.iLDeN)) THEN
!
!  The POM flux reaching the seafloor is partially resuspended in the
!  lower water column level within the small detritus compartment,
!  assuming fragmentation of large particles.
!  Druon 1999 (pg 32) assumes the resuspension rate (ReSuspR:
!  nondimensional) is a function of the bottom friction velocity and
!  the critical friction velocity for resuspension of Particulate
!  Organic Matter: ReSuspR = (Ustarb/Ustar_crit)^2 where Ustarb is
!  computed from the model bottom stress, bustr. The critical stress
!  proposed by Peterson 1999 (Aquacultural Engineering, 21, 85-111)
!  is 0.01 Pa equals Ustar_crit=3.1E-3 m s-1 for seawater rho = 1027.
!  The factor 9.61E-6 is 3.1E-3^2.
!
!  The sinking flux FC(i,0) is in units millimole m-2 because it is
!  concentration times layer thickness and has been integrated over
!  model dt. When converting this to a quantity to be added to a
!  tracer variable (concentration in millimole m-3) must divide
!  by layer thickness.
!
                Ustarb=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+ &
     &                 (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
                ReSuspR=MIN(Ustarb**2/9.61E-6_r8, 1.0_r8) ! Druon et al.
!                                                           2010
!  The particulate flux reaching the seafloor that will
!  be remineralized/denitrified.
!
                cff1=FC(i,0)*Hz_inv(i,1)*(1.0_r8-ReSuspR) ! mmol m-3
!
!  The complimentary resuspended fraction goes back to small detritus.
!
                cff5=FC(i,0)*Hz_inv(i,1)*ReSuspR
!
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  Bio(i,1,iSDeN)=Bio(i,1,iSDeN)+cff5
                END IF
                IF ((ibio.eq.iSDeC).or.(ibio.eq.iLDeC)) THEN
                  Bio(i,1,iSDeC)=Bio(i,1,iSDeC)+cff5
                END IF
                IF (ibio.eq.iPhyt) THEN
                  Bio(i,1,iSDeC)=Bio(i,1,iSDeC)+cff5*PhyCN(ng)
                END IF
              END IF
!
!  Nitrogen burial efficiency.
!
              N_BURIAL : IF ((ibio.eq.iPhyt).or.                        &
     &                       (ibio.eq.iSDeN).or.                        &
     &                       (ibio.eq.iLDeN)) THEN
!
!  Burial efficiency (Cbe) is computed with a parameterization based on
!  Henrichs and Reeburgh (1987), "Anerobic mineralization of marine
!  sediment organic matter: Rates and the role of anerobic processes in
!  the oceanic carbon economy", Geomicrobiology Journal, 5, 191-237.
!  Combining their equations (2) and (3) gives a relation
!  for the fractional burial efficiency, E = 0.023*F^0.5797, as a
!  function of organic carbon flux, F (given in units gC m-2 yr-1).
!
!  To estimate particulate nitrogen burial efficiency, we convert the
!  particulate nitrogen fluxes to equivalent carbon fluxes by assuming
!  a fixed CN ratio, CNbur (Gelinas et al. 2001), limited by a
!  maximum value.
!
!  Note that Henrichs empirical relation is for total particulate
!  flux combining all detritus and phytoplankton, whereas here we
!  consider each constituent separately in the loop over indx. The
!  Henrichs relation is nonlinear, so this approach incurs an error.
!
!  The sinking flux FC(i,0) is in units millimole m-2 because it is
!  concentration times layer thickness and has been integrated over
!  model dt. cff1 is FC(i,0)/Hz so carbon flux in gC m-2 year-1 (fac1)
!  is computed as follows.
!
!  JW: In original code the same CNbur of 9.3 was used for phytoplankton
!  which gives inconsistent burial efficiency for phytoplankton C and N.
!
                IF (ibio.eq.iPhyt) THEN
                  CNbur=PhyCN(ng)
                ELSE
                  CNbur=9.3_r8
                ENDIF
                cff5=0.01_r8           ! 1 percent of flux goes to DON
                fac1=CNbur*                                             &
     &               12.0_r8/1000.0_r8*                                 &
     &               cff1*Hz(i,j,1)*(365.0_r8/dtdays) ! gC m-2 year-1
                Cbe=MIN(0.50_r8,0.023_r8*fac1**0.5797_r8 )
#ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iNbur)=DiaBio2d(i,j,iNbur)+                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              Cbe*cff1*Hz(i,j,1)*fiter
#endif
!
!  cff1 will now be the unburied flux that is to be remineralized
!  Must multiply cff1 by Hz to convert to flux in meter-2 when saving
!  the denitrification diagnostic term. Leave in meter-3 when returning
!  C or N to dissolved tracer concentrations.
!
                cff1=(1.0_r8-Cbe)*cff1                        ! mmol m-3
!
                TS=LOG((298.15_r8-Bio(i,1,itemp))/                      &
     &             (273.15_r8+Bio(i,1,itemp)))

                AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+      &
     &             Bio(i,1,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,1,isalt)*Bio(i,1,isalt)
                O2satB=l2mol*EXP(AA)
                LBO2=(K_BO2*(O2satB-Bio(i,1,iOxyg)))/                   &
     &                  (O2satB*(Bio(i,1,iOxyg)+K_BO2))

#ifdef NET_SULFATE_REDUCTION
!
!  Following Da et al. (2024).
!
!  Assume 30% of OM is oxidized by O2 and NO3, but 70% of OM is oxidized
!  by sulfate (Mackin and Swider, 1989). Note: cff2=4/16, cff5=1/100.
!
                fac=cff1*(1.0_r8+3.0_r8*LBO2)*0.30_r8
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff2*fac
                Bio(i,1,iTAlk)=Bio(i,1,iTAlk)+cff2*fac
                Bio(i,1,iDON_)=Bio(i,1,iDON_)+cff5*fac
!
!  70% of OM is oxidized by sulfate, consume NH4 and NO3.
!  NH4 production at the sediment (16NH3 part) minus NH4 consumption
!  at the sediment boundary (net sulfate reduction part).
!
                fac=106.0_r8/96.0_r8*(1.0_r8-LBO2)
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+                          &
     &                         cff1*0.70_r8*(1.0_r8-fac)
                Bio(i,1,iNH4_)=MAX(0.0_r8, Bio(i,1,iNH4_))
!
!  NO3 consumption at the sediment boundary.
!
                Bio(i,1,iNO3_)=Bio(i,1,iNO3_)-                          &
     &                         cff1*0.70_r8*                            &
     &                         106.0_r8/160.0_r8*(1.0_r8-LBO2)
                Bio(i,1,iNO3_)=MAX(0.0_r8, Bio(i,1,iNO3_))
#else
                fac=cff1*(1.0_r8+3.0_r8*LBO2)
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff2*fac
                Bio(i,1,iTAlk)=Bio(i,1,iTAlk)+cff2*fac
                Bio(i,1,iDON_)=Bio(i,1,iDON_)+cff5*fac
#endif
#ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iDNIT)=DiaBio2d(i,j,iDNIT)+                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              (1.0_r8-(cff2+cff5)*                &
     &                                      (1.0_r8+3.0_r8*LBO2))*      &
# ifdef NET_SULFATE_REDUCTION
     &                              0.30_r8*                            &
# endif
     &                              cff1*Hz(i,j,1)*fiter
#endif
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-                          &
     &                         cff1*cff3*                               &
#ifdef NET_SULFATE_REDUCTION
     &                         0.30_r8*                                 &
#endif
     &                         (1.0_r8-LBO2)
#ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iSoxy)=DiaBio2d(i,j,iSoxy)+                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              cff1*cff3*(1.0_r8-LBO2)*            &
# ifdef NET_SULFATE_REDUCTION
     &                              0.30_r8*                            &
# endif
     &                              Hz(i,j,1)*fiter
#endif
              END IF N_BURIAL
!
!  Carbon burial efficiency.
!
              C_BURIAL : IF ((ibio.eq.iPhyt).or.                        &
     &                       (ibio.eq.iSDeC).or.                        &
     &                       (ibio.eq.iLDeC)) THEN
!
!  Here, cff1 currently represents the non-resuspended, non-buried
!  phytoplankton flux. However, the Cbe calculation below requires
!  the non-resuspended phytoplankton flux. Therefore, cff1 needs to
!  be rescaled:
!
                IF (ibio.eq.iPhyt) THEN
                  cff1=cff1/(1.0_r8-Cbe)
                END IF
!
!  In the iPhyt case we are repeating the Cbe calculation above.
!  CNbur coefficient not required here because we are already in carbon
!  units except for phytoplankton, in which case use PhyCN ratio.
!
                fac1=12.0_r8/1000.0_r8*                                 &
     &               cff1*Hz(i,j,1)*(365.0_r8/dtdays) ! gC m-2 year-1
                IF (ibio.eq.iPhyt) THEN
                  fac1=fac1*PhyCN(ng)
                ENDIF
                Cbe=MIN(0.50_r8,0.023_r8*fac1**0.5797_r8 )
!
!  Factor for converting phytoplankton nitrogen flux to carbon.
!
                IF (ibio.eq.iPhyt) THEN
                  fac2=PhyCN(ng)
                ELSE
                  fac2=1.0_r8
                ENDIF

#ifdef DIAGNOSTICS_BIO
!
!  The flux of carbon that reaches the seafloor.
!
                DiaBio2d(i,j,iCbot)=DiaBio2d(i,j,iCbot)+                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              fac2*cff1*Hz(i,j,1)*fiter
!
!  The flux of carbon that is buried.
!
                DiaBio2d(i,j,iCbur)=DiaBio2d(i,j,iCbur)+                &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
# endif
     &                              Cbe*fac2*cff1*Hz(i,j,1)*fiter
#endif
!
!  The fraction of the flux that is not buried will be remineralized.
!

                cff1=(1.0_r8-Cbe)*cff1                        ! mmol m-3
!
                IF (ibio.eq.iPhyt) THEN
                  fac2=PhyCN(ng)    ! phytoplankton nitrogen to carbon
                ELSE
                  fac2=1.0_r8
                ENDIF
                cff5=0.01_r8        ! 1% of flux goes to DOC (as for DON)

#ifdef NET_SULFATE_REDUCTION
!
!  Following Da et al. (2024).
!
                Bio(i,1,iDOC_)=Bio(i,1,iDOC_)+                          &
     &                         fac2*cff1*cff5*0.30_r8
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+                          &
     &                         fac2*cff1*(1.0_r8-cff5)*0.30_r8
!
!  70% of the organic matter will be remineralized by net sulfate
!  reduction. Here, we use organic carbon for all calculations. This
!  net sulfate reduction reaction consumes NH4 and NO3 and O2. A
!  replacement might be needed for LBO2 because KBO2 used for LBO2
!  calculation is for bottom denitrification switch.
!
                Bio(i,1,iTAlk)=Bio(i,1,iTAlk)+                          &
     &                         fac2*cff1*0.70_r8*                       &
     &                         (106.0_r8*(1.0_r8-1.1_r8*2.0_r8/         &
     &                          3.0_r8*(1.0_r8-LBO2))+15.0_r8)/         &
     &                         106.0_r8
                fac=fac2*cff1*0.70_r8
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+                          &
     &                         fac
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-                          &
     &                         fac*(2.0_r8/3.0_r8*(1.0_r8-LBO2))
!
!  Add additional O2 consumption via sulfate reduction.
!
# ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iSoxy)=DiaBio2d(i,j,iSoxy)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              fac*(2.0_r8/3.0_r8*(1.0_r8-LBO2))*  &
     &                              Hz(i,j,1)*fiter
# endif
#else
                fac=fac2*cff1
                Bio(i,1,iDOC_)=Bio(i,1,iDOC_)+fac*cff5
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+fac*(1.0_r8-cff5)
#endif
              END IF C_BURIAL
            END DO

          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012).
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iTIC_)=MIN(Bio(i,k,iTIC_), 3000.0_r8)
            Bio(i,k,iTIC_)=MAX(Bio(i,k,iTIC_), 1.0_r8)
          END DO
        END DO
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
#ifdef MASKING
              cff=cff*rmask(i,j)
# ifdef WET_DRY
              cff=cff*rmask_wet(i,j)
# endif
#endif
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP
!
      RETURN
      END SUBROUTINE ecb_tile

#ifdef pCO2_RZ
      SUBROUTINE pCO2_water_RZ (Istr, Iend,                             &
     &                          LBi, UBi, LBj, UBj, IminS, ImaxS,       &
     &                          j, DoNewton,                            &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          T, S, TIC, TAlk, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4b=0,           !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2= ppmv  (DoNewton=0)                              !
!                pCO2= ppmv  (DoNewton=1)                              !
!                                                                      !
!  This subroutine was adapted by Katja Fennel (Nov 2005) from         !
!  Zeebe and Wolf-Gladrow (2001).                                      !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Zeebe, R.E. and D. Wolf-Gladrow,  2005:  CO2 in Seawater:         !
!      Equilibrium, kinetics, isotopes, Elsevier Oceanographic         !
!      Series, 65, pp 346.                                             !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif

      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: scl, sqrtS
      real(r8) :: borate, alk, dic
      real(r8) :: ff, K1, K2, K12, Kb, Kw
      real(r8) :: p5, p4, p3, p2, p1, p0
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
# ifdef pCO2_RZ_MILLERO_2010
      real(r8) :: pkmi, a_mi, b_mi, c_mi
# endif
# ifdef pCO2_RZ_CAIWANG_1998
      real(r8) :: pkcw, fhcw
# endif
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
# ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
# endif
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        scl=S(i)/1.80655_r8

        alk= TAlk(i)*0.000001_r8
        dic = TIC(i)*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carbonic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
# if defined pCO2_RZ_MILLERO_2010

!  Use K1,K2 of Millero, 2010, Marine and Freshwater Research, v.61,
!  p.139-142. The coefficients are on the seawater scale (SWS).
!-----------------------------------------------------------------------
!
        a_mi=13.40380_r8*sqrtS+                                         &
     &       S(i)*(0.032060_r8-5.242e-5_r8*S(i))
        b_mi=-530.65900_r8*sqrtS-5.821000_r8*S(i)
        c_mi=-2.06640_r8*sqrtS
        pkmi=-126.34048_r8+6320.813_r8*invTk+19.568224_r8*logTk
        pkmi=pkmi+a_mi+b_mi*invTk+c_mi*logTk
        K1=10.0_r8**(-pkmi)
!
        a_mi=21.37280_r8*sqrtS+                                         &
     &       S(i)*(0.121800_r8-3.688e-4_r8*S(i))
        b_mi=-788.28900_r8*sqrtS-19.189000_r8*S(i)
        c_mi=-3.37400_r8*sqrtS
        pkmi=-90.18333_r8+5143.692_r8*invTk+14.613358_r8*logTk
        pkmi=pkmi+a_mi+b_mi*invTk+c_mi*logTk
        K2=10.0_r8**(-pkmi)

# elif defined pCO2_RZ_CAIWANG_1998

!  Use K1,K2 of Cai and Wang, 1998, Limnology and Oceanography, p. 661.
!  The coefficients are on NBS scale, need to use "fhcw" to convert to
!  SWS. "fhcw" calculation is from CO2SYS.m (Takahashi et al., 1982).
!-----------------------------------------------------------------------
!
        fhcw=1.2948_r8-2.036e-3_r8*Tk+                                  &
     &       (4.607e-4_r8-1.475e-6_r8*Tk)*S(i)*S(i)
        pkcw=200.1_r8  *invTk+0.322_r8
        pkcw=3404.71_r8*invTk+0.032786_r8*Tk-14.8435_r8-                &
     &       7.1692e-2_r8*pkcw*sqrtS+2.1487e-3_r8*S(i)
        K1=10.0_r8**(- pkcw)     ! This is on NBS scale
        K1=K1/fhcw               ! Convert to SWS scale
!
        pkcw=-129.24_r8*invTk+1.4381_r8
        pkcw=2902.39_r8*invTk+0.02379_r8*Tk-6.498_r8-                   &
     &       0.3191_r8*pkcw*sqrtS+0.0198_r8*S(i)
        K2=10.0_r8**(-pkcw)      ! This is on NBS scale
        K2=K2/fhcw               ! Convert to SWS scale
# else

!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
# endif
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Solve for h in fifth-order polynomial. First calculate
!  polynomial coefficients.
!
        K12 = K1*K2

        p5 = -1.0_r8
        p4 = -alk-Kb-K1
        p3 = dic*K1-alk*(Kb+K1)+Kb*borate+Kw-Kb*K1-K12
        p2 = dic*(Kb*K1+2*K12)-alk*(Kb*K1+K12)+Kb*borate*K1             &
     &       +(Kw*Kb+Kw*K1-Kb*K12)
        p1 = 2.0_r8*dic*Kb*K12-alk*Kb*K12+Kb*borate*K12                 &
     &       +Kw*Kb*K1+Kw*K12
        p0 = Kw*Kb*K12
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
!
          DO Inewton=1,InewtonMax
!
!  Evaluate f([H+]) = p5*x^5+...+p1*x+p0
!
            fn=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
!
!  Evaluate derivative, df([H+])/dx:
!
!     df= d(fn)/d(X)
!
            df=(((5*p5*X+4*p4)*X+3*p3)*X+2*p2)*X+p1
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0) THEN
               EXIT BRACK_IT
            ELSE
               ftest=fni(1)/fni(3)
               IF (ftest.gt.0) THEN
                 X_hi=X_mid
               ELSE
                 X_lo=X_mid
               END IF
               X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

# ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
# endif

      END DO I_LOOP
!
      RETURN
      END SUBROUTINE pCO2_water_RZ
#else
      SUBROUTINE pCO2_water (Istr, Iend,                                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j, DoNewton,                 &
# ifdef MASKING
     &                       rmask,                                     &
# endif
     &                       T, S, TIC, TAlk, PO4b, SiO3, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     PO4b       Inorganic phosphate (millimol/m3).                    !
!     SiO3       Inorganic silicate (millimol/m3).                     !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4b=0,           !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2=0.35074945E+03 ppmv  (DoNewton=0)                !
!                pCO2=0.35073560E+03 ppmv  (DoNewton=1)                !
!                                                                      !
!  This subroutine was adapted by Mick Follows (Oct 1999) from OCMIP2  !
!  code CO2CALC. Modified for ROMS by Hernan Arango (Nov 2003).        !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: PO4b
      real(r8), intent(in) :: SiO3

      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: SO4, scl, sqrtS, sqrtSO4
      real(r8) :: alk, dic, phos, sili
      real(r8) :: borate, sulfate, fluoride
      real(r8) :: ff, K1, K2, K1p, K2p, K3p, Kb, Kf, Ks, Ksi, Kw
      real(r8) :: K12, K12p, K123p, invKb, invKs, invKsi
      real(r8) :: A, A2, B, B2, C, dA, dB
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
# ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
# endif
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        SO4=19.924_r8*S(i)/(1000.0_r8-1.005_r8*S(i))
        sqrtSO4=SQRT(SO4)
        scl=S(i)/1.80655_r8

        alk=TAlk(i)*0.000001_r8
        dic=TIC(i)*0.000001_r8
        phos=PO4b*0.000001_r8
        sili=SiO3*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute first (K1p), second (K2p), and third (K3p) dissociation
!  constant of phosphoric acid:
!
!           K1p = [H][H2PO4]/[H3PO4]
!           K2p = [H][HPO4]/[H2PO4]
!           K3p = [H][PO4]/[HPO4]
!
!  From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
!  With footnote using data from Millero (1974).
!-----------------------------------------------------------------------
!
        K1p=EXP(115.525_r8-                                             &
     &          invTk*4576.752_r8-                                      &
     &          logTk*18.453_r8+                                        &
     &          sqrtS*(0.69171_r8-invTk*106.736_r8)-                    &
     &          S(i)*(0.01844_r8+invTk*0.65643_r8))
        K2p=EXP(172.0883_r8-                                            &
     &          invTk*8814.715_r8-                                      &
     &          logTk*27.927_r8+                                        &
     &          sqrtS*(1.3566_r8-invTk*160.340_r8)-                     &
     &          S(i)*(0.05778_r8-invTk*0.37335_r8))
        K3p=EXP(-18.141_r8-                                             &
     &          invTk*3070.75_r8+                                       &
     &          sqrtS*(2.81197_r8+invTk*17.27039_r8)-                   &
     &          S(i)*(0.09984_r8+invTk*44.99486_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
!  From Millero (1995; page 671) using data from Yao and Millero (1995).
!-----------------------------------------------------------------------
!
        Ksi=EXP(117.385_r8-                                             &
     &          invTk*8904.2_r8-                                        &
     &          logTk*19.334_r8+                                        &
     &          sqrtSO4*(3.5913_r8-invTk*458.79_r8)-                    &
     &          SO4*(1.5998_r8-invTk*188.74_r8-                         &
     &               SO4*(0.07871_r8-invTk*12.1652_r8))+                &
     &          LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!------------------------------------------------------------------------
!  Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
!  From Dickson (1990, J. chem. Thermodynamics 22, 113)
!------------------------------------------------------------------------
!
        Ks=EXP(141.328_r8-                                              &
     &         invTk*4276.1_r8-                                         &
     &         logTk*23.093_r8+                                         &
     &         sqrtSO4*(324.57_r8-invTk*13856.0_r8-logTk*47.986_r8-     &
     &                  SO4*invTk*2698.0_r8)-                           &
     &         SO4*(771.54_r8-invTk*35474.0_r8-logTk*114.723_r8-        &
     &              SO4*invTk*1776.0_r8)+                               &
     &         LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
!  From Dickson and Riley (1979) -- change pH scale to total.
!-----------------------------------------------------------------------
!
        Kf=EXP(-12.641_r8+                                              &
     &         invTk*1590.2_r8+                                         &
     &         sqrtSO4*1.525_r8+                                        &
     &         LOG(1.0_r8-0.001005_r8*S(i))+                            &
     &         LOG(1.0_r8+0.1400_r8*scl/(96.062_r8*Ks)))
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974), sulfate (Morris
! and Riley, 1966), and fluoride (Riley, 1965).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
        sulfate=0.14_r8*scl/96.062_r8
        fluoride=0.000067_r8*scl/18.9984_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          DO Inewton=1,InewtonMax
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
            X2=X*X
            X3=X2*X
            invX=1.0_r8/X
            invX2=1.0_r8/X2

            A=X*(K12p+X*(K1p+X))
            B=X*(K1+X)+K12
            C=1.0_r8/(1.0_r8+sulfate*invKs)

            A2=A*A
            B2=B*B
            dA=X*(2.0_r8*K1p+3.0_r8*X)+K12p
            dB=2.0_r8*X+K1
!
!  Evaluate f([H+]):
!
!     fn=HCO3+CO3+borate+OH+HPO4+2*PO4+H3PO4+silicate+Hfree+HSO4+HF-TALK
!
            fn=dic*K1*(X+2.0_r8*K2)/B+                                  &
     &         borate/(1.0_r8+X*invKb)+                                 &
     &         Kw*invX+                                                 &
     &         phos*(K12p*X+2.0_r8*K123p-X3)/A+                         &
     &         sili/(1.0_r8+X*invKsi)-                                  &
     &         X*C-                                                     &
     &         sulfate/(1.0_r8+Ks*invX*C)-                              &
     &         fluoride/(1.0_r8+Kf*invX)-                               &
     &         alk
!
!  Evaluate derivative, f(prime)([H+]):
!
!     df= d(fn)/d(X)
!
            df=dic*K1*(B-dB*(X+2.0_r8*K2))/B2-                          &
     &         borate/(invKb*(1.0+X*invKb)**2)-                         &
     &         Kw*invX2+                                                &
     &         phos*(A*(K12p-3.0_r8*X2)-dA*(K12p*X+2.0_r8*K123p-X3))/A2-&
     &         sili/(invKsi*(1.0_r8+X*invKsi)**2)+                      &
     &         C+                                                       &
     &         sulfate*Ks*C*invX2/((1.0_r8+Ks*invX*C)**2)+              &
     &         fluoride*Kf*invX2/((1.0_r8+Kf*invX)**2)
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
              X2=X*X
              X3=X2*X
              invX=1.0_r8/X

              A=X*(K12p+X*(K1p+X))+K123p
              B=X*(K1+X)+K12
              C=1.0_r8/(1.0_r8+sulfate*invKs)

              A2=A*A
              B2=B*B
              dA=X*(K1p*2.0_r8+3.0_r8*X2)+K12p
              dB=2.0_r8*X+K1
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=dic*(K1*X+2.0_r8*K12)/B+                       &
     &                   borate/(1.0_r8+X*invKb)+                       &
     &                   Kw*invX+                                       &
     &                   phos*(K12p*X+2.0_r8*K123p-X3)/A+               &
     &                   sili/(1.0_r8+X*invKsi)-                        &
     &                   X*C-                                           &
     &                   sulfate/(1.0_r8+Ks*invX*C)-                    &
     &                   fluoride/(1.0_r8+Kf*invX)-                     &
     &                   alk
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0.0_r8) THEN
              EXIT BRACK_IT
            ELSE
              ftest=fni(1)/fni(3)
              IF (ftest.gt.0.0) THEN
                X_hi=X_mid
              ELSE
                X_lo=X_mid
              END IF
              X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

# ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
# endif

      END DO I_LOOP
!
      RETURN
      END SUBROUTINE pCO2_water
#endif
      END MODULE biology_mod
