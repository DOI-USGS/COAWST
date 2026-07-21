      MODULE ice_thermo_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group           Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!================================================== Hernan G. Arango ===
!                                                                      !
!  It computes the ice thermodynamic growth and decay term based on    !
!  Mellor and Kantha (1989) and Parkinson and Washington (1979).       !
!                                                                      !
!  It computes heat fluxes and ice production rates:                   !
!                                                                      !
!    Fi(i,j,icW_ai) = -(qai(i,j) - qi2(i,j)) /(hfus1(i,j)*RhoSW)       !
!                                                                      !
!  and updates the internal ice temperature:                           !
!                                                                      !
!    Si(i,j,linew,isTice) = Si(i,j,linew,isTice) + dtice*RHS           !
!                                                                      !
!  Required atmospheric fields:                                        !
!                                                                      !
!    Qnet_ai(:,:)          net heat flux at the air-ice interface      !
!    Qnet_ao(:,:)          net heat flux at the air-ocean interface    !
!    rain(:,:)             rain fall rate                              !
!    snow(:,:)             snow fall rate                              !
!    sustr(:,:)            surface u-wind stress                       !
!    svstr(:,:)            surface v-wind stress                       !
!                                                                      !
!  Required ocean fields:                                              !
!                                                                      !
!    dtice                 ice kernel timestep                         !
!    rmask(:,:)            land/sea mask                               !
!    t(:,:,N(ng),:,:)      surface ocean temperature and salinity      !
!                                                                      !
!  Required ice variables:                                             !
!                                                                      !
!    Fi(:,:,icQcon)        gradient coefficient for heat conductivity  !
!    Fi(:,:,icQrhs)        downward heat conductivity term             !
!    Fi(:,:,icS0mk)        salinity of molecular sublayer under ice    !
!    Fi(:,:,icT0mk)        temperature of molecular sublayer under ice !
!    Fi(:,:,icIsst)        temperature at snow-air interface           !
!    Fi(:,:,icW_ai)        rate of melt/freeze at air-ice interface    !
!    Fi(:,:,icW_ao)        rate of melt/freeze at air-ocena interface  !
!    Fi(:,:,icW_fr)        rate of ice accretion by frazil growth      !
!    Fi(:,:,icW_io)        rate of melt/freeze at ice-ocean interface  !
!    Fi(:,:,icW_ro)        rate of melt/freeze runoff into ocean       !
!    Fi(:,:,icIOmf)        ice-ocean mass flux                         !
!                                                                      !
!    Si(:,:,linew,isAice)  ice concentration                           !
!    Si(:,:,linew,isIage)  ice age                                     !
!    Si(:,:,linew,isEnth)  scaled perturbation ice heat content        !
!    Si(:,:,linew,isHice)  ice thickness, ice mass (divided by area)   !
!    Si(:,:,linew,isHsno)  snow thickness, mass snow per area          !
!    Si(:,:,linew,isHmel)  melt water thickness on ice                 !
!    Si(:,:,linew,isTice)  ice interior temperature (ice layer middle) !
!                                                                      !
!  Relevant Internal variables:                                        !
!                                                                      !
!    brnfr(:,:)            brine fraction                              !
!    hfus1(:,:)            latent heat of fusion (L_o or L_3)          !
!    ice_thick(:,:)        ice thickness                               !
!    qai(:,:)              upward air-ice heat flux                    !
!    qio(:,:)              upward ice-ocean heat flux                  !
!    qi2(:,:)              heat flux in ice                            !
!    sice(:,:)             ice salinity                                !
!    snow_thick(:,:)       snow thickness                              !
!    t2(:,:)               temperature at ice/snow interface           !
!    wsm(:,:)              snow melting rate                           !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Mellor, G.L. and L. Kantha, 1989: An Ice-Ocean Coupled Model,     !
!      J. Geophys. Res., 94, 10937-10954.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
#ifdef AICLM_NUDGING
      USE mod_clima
#endif
#ifdef ICE_SHOREFAST
      USE mod_coupling
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod,       ONLY : bc_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
      PUBLIC  :: ice_thermo
      PRIVATE :: ice_thermo_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_thermo (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
#endif
      CALL ice_thermo_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), liold(ng), linew(ng),             &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
#ifdef ICE_SHOREFAST
     &                      GRID(ng) % h,                               &
     &                      COUPLING(ng) % Zt_avg1,                     &
#endif
#ifdef AICLM_NUDGING
     &                      CLIMA(ng) % aiclm,                          &
     &                      CLIMA(ng) % hiclm,                          &
     &                      CLIMA(ng) % AInudgcof,                      &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % Qnet_ai,                       &
     &                      FORCES(ng) % Qnet_ao,                       &
     &                      FORCES(ng) % snow,                          &
     &                      FORCES(ng) % rain,                          &
     &                      FORCES(ng) % stflx,                         &
     &                      ICE(ng) % Fi,                               &
     &                      ICE(ng) % Si)
#ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE ice_thermo
!
!***********************************************************************
      SUBROUTINE ice_thermo_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, liold, linew,                   &
#ifdef MASKING
     &                            rmask,                                &
#endif
#ifdef WET_DRY
     &                            rmask_wet,                            &
#endif
#ifdef ICESHELF
     &                            zice,                                 &
#endif
     &                            z_r, z_w,                             &
#ifdef ICE_SHOREFAST
     &                            h, Zt_avg1,                           &
#endif
#ifdef AICLM_NUDGING
     &                            aiclm, hiclm, AInudgcof,              &
#endif
     &                            t,                                    &
     &                            sustr, svstr,                         &
     &                            Qnet_ai, Qnet_ao,                     &
     &                            snow,                                 &
     &                            rain,                                 &
     &                            stflx,                                &
     &                            Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, liold, linew
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
# endif
# ifdef AICLM_NUDGING
      real(r8), intent(in) :: aiclm(LBi:,LBj:)
      real(r8), intent(in) :: hiclm(LBi:,LBj:)
      real(r8), intent(in) :: AInudgcof(LBi:,LBj:)
# endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: Qnet_ai(LBi:,LBj:)
      real(r8), intent(in) :: Qnet_ao(LBi:,LBj:)
      real(r8), intent(in) :: snow(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
# endif
# ifdef AICLM_NUDGING
      real(r8), intent(in) :: aiclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hiclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: AInudgcof(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Qnet_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Qnet_ao(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: snow(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
#endif

! Local variable definitions
!
      logical :: IceCavity
!
      integer :: i, j
!
      real(r8), parameter :: AlphIc = 2.034_r8          ! [W m-1 K-1]
      real(r8), parameter :: AlphSn = 0.31_r8           ! [W m-1 K-1]
      real(r8), parameter :: Cp_i = 2093.0_r8           ! [J kg-1 K-1]
      real(r8), parameter :: Cp_w = 3990.0_r8           ! [J kg-1 K-1]
      real(r8), parameter :: eps = 1.0E-4_r8            ! zero division
      real(r8), parameter :: frln = -0.0543_r8          ! [psu C-1]
      real(r8), parameter :: hfus = 3.347E+5_r8         ! [J kg-1]
      real(r8), parameter :: kappa = 0.4_r8             ! von Karman
      real(r8), parameter :: nu = 1.8E-6_r8             ! m2/s
      real(r8), parameter :: prs = 2432.0_r8            ! S Schmidt Num.
      real(r8), parameter :: prt = 13.0_r8              ! T Prandtl Num.
      real(r8), parameter :: RhoCpr = 0.2442754E-6_r8   ! [m s2 K kg-1]
      real(r8), parameter :: RhoSW = 1026.0_r8          ! [kg m-3]
      real(r8), parameter :: sice_ref = 3.2_r8          ! [psu]
      real(r8), parameter :: tpr = 0.85_r8              ! Turb. Prandtl
      real(r8), parameter :: ykf = 3.14                 ! Yaglom/Kader
      real(r8), parameter :: z0ii = 0.02_r8             ! ice roughness
!
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: d1, d2i, d3, dztop, fac_shflx
      real(r8) :: ai_tmp, corfac, cot, delta_mi
      real(r8) :: hicehinv, hstar, mi_old, phi
      real(r8) :: Qsur, rno, termt, terms, tfrz, tfz
      real(r8) :: xwai, xtot, z0, zdz0, xmelt
#ifdef ICE_SHOREFAST
      real(r8) :: clear, fac_sf, hh
#endif
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: alph
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: brnfr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: b2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: cht
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Coa
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: hfus1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qai
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qio
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qi2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: salt_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: temp_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: t2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utau
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ws
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wsm

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initalize.
!-----------------------------------------------------------------------
!
!  Set bulk sensible transfer coefficient (d1), latent heat transfer
!  over ice (d2i), and Stefan-Boltzan constant times surface
!  emissivity (d3).
!
      d1 =AirRho(ng)*spec_heat_air*trans_coeff            ! [J/(K m3)]
      d2i=AirRho(ng)*sublimation*trans_coeff              ! [J/(K m3)]
      d3 =StefBo*ice_emiss                                ! [W/(K^4 m2)]
!
!  Extract ocean surface temperature and salinity. Compute surface
!  level thickness.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          temp_top(i,j)=t(i,j,N(ng),nrhs,itemp)
          salt_top(i,j)=t(i,j,N(ng),nrhs,isalt)
          salt_top(i,j)=MIN(MAX(0.0_r8, salt_top(i,j)), 40.0_r8)
        END DO
      END DO
!
!  Compute squared-root of surface wind stress magnitude.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          utau(i,j)=SQRT(SQRT((0.5_r8*(sustr(i  ,j)+                    &
                                       sustr(i+1,j)))**2+               &
     &                        (0.5_r8*(svstr(i,j  )+                    &
     &                                 svstr(i,j+1)))**2))
          utau(i,j)=MAX(utau(i,j), 1.0E-4_r8)
        END DO
      END DO
!
!  Compute snow thickness, ice thickness, brine fraction, and
!  thermal conductivity.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sice(i,j)=MIN(sice_ref, salt_top(i,j))
          ice_thick(i,j)=0.05_r8+                                       &
     &                   Si(i,j,linew,isHice)/                          &
     &                   (Si(i,j,linew,isAice)+eps)
          snow_thick(i,j)=Si(i,j,linew,isHsno)/                         &
                          (Si(i,j,linew,isAice)+eps)
          brnfr(i,j)=frln*sice(i,j)/(Si(i,j,linew,isTice)-eps)
          brnfr(i,j)=MIN(brnfr(i,j),0.2_r8)
          brnfr(i,j)=MAX(brnfr(i,j),0.0_r8)
          alph(i,j)=AlphIc*MAX(1.0_r8-1.2_r8*brnfr(i,j), 0.25_r8)
          cff=(Si(i,j,linew,isHice)/1.0_r8)**2
          corfac=1.0_r8/(0.5_r8*(1.0_r8+EXP(-cff)))
          alph(i,j)=alph(i,j)*corfac
          Coa(i,j)=2.0_r8*alph(i,j)*snow_thick(i,j)/                    &
     &             (AlphSn*ice_thick(i,j))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Solve for temperature at the top of the ice layer.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Gradient coefficient for heat conductivity term.
!
          b2d(i,j)=2.0_r8*alph(i,j)/(ice_thick(i,j)*(1.0_r8+Coa(i,j)))
          Fi(i,j,icQcon)=Fi(i,j,icQcon)+                                &
     &                   b2d(i,j)
!
!  Downward conductivity term, assuming the ocean at the freezing point
!  (convert ice temperature to Kelvin).
!
          IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
            Fi(i,j,icQrhs)=Fi(i,j,icQrhs)+                              &
     &                     b2d(i,j)*(Si(i,j,linew,isTice)+273.15_r8)
!
!  Compute temperature at the snow/ice interface (convert to Celsius).
!
            Fi(i,j,icIsst)=(Fi(i,j,icQrhs)/Fi(i,j,icQcon))-273.15_r8
!
!  Bound value at zero Celsius for stability. It can occassionably be
!  unstable and take "icIsst" in the wrong direction.
!
            Fi(i,j,icIsst)=MIN(MAX(Fi(i,j,icIsst),-45.0_r8), 0.0_r8)
          ELSE
            Fi(i,j,icIsst)=temp_top(i,j)
          END IF
        END DO
      END DO
!
!  Calculate new interior ice temperature.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  The "cot" calculation derives from an assumption of a linear
!  relationship between ice temperature and salinity.  As "isTice"
!  approaches zero Celsius, "cot" goes to infinity and any change
!  in "isTice" in the timestepping below is stymied, so keep the
!  denominator here below zero (SMD).
!
          IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
            cot=-frln*sice(i,j)*hfus/                                   &
     &          (MIN(Si(i,j,linew,isTice), frln*sice_ref))**2+Cp_i
            cff1=IceRho(ng)*cot*ice_thick(i,j)**2
            cff2=Fi(i,j,icIsst)-(2.0_r8+Coa(i,j))*Si(i,j,linew,isTice)
            cff3=1.0_r8+Coa(i,j)
            Si(i,j,linew,isTice)=Si(i,j,linew,isTice)+                  &
     &                           dtice(ng)*                             &
     &                           (2.0_r8*alph(i,j)/cff1*                &
     &                            (Fi(i,j,icT0mk)+cff2/cff3))
            Si(i,j,linew,isTice)=MAX(Si(i,j,linew,isTice), -35.0_r8)
!
!  Ensure that "isTice" remains below "frln*sice_ref" whenever ice is
!  present. Also ensure that it remains below the maximum of either
!  the ice surface temperature or the water temperature below. This
!  is imperfect as warmer ice can advect from elsewhere and alter
!  its heat content, or thick ice might actually have a warmer
!  interior temperature than either the surface or bottom. But in
!  general it does not make sense that ice forms at a warmer
!  temperature than exists in any of its surroundings (SMD).
!
            Si(i,j,linew,isTice)=MIN(Si(i,j,linew,isTice),              &
     &                               frln*sice_ref)
            Si(i,j,linew,isTice)=MIN(Si(i,j,linew,isTice),              &
     &                               MAX(Fi(i,j,icIsst),                &
     &                                   Fi(i,j,icT0mk)))
          ELSE
            Si(i,j,linew,isTice)=temp_top(i,j)
          END IF
        END DO
      END DO
!
!  Calculate associated heat fluxes.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
            hicehinv=1.0_r8/(0.5_r8*ice_thick(i,j))
            t2(i,j) =(Fi(i,j,icIsst)+Coa(i,j)*Si(i,j,linew,isTice))/    &
     &               (1.0_r8+Coa(i,j))
            Qi2(i,j)=alph(i,j)*                                         &
     &               (Si(i,j,linew,isTice)-t2(i,j))*hicehinv
            Qio(i,j)=alph(i,j)*                                         &
     &               (Fi(i,j,icT0mk)-Si(i,j,linew,isTice))*hicehinv
          END IF
          Qai(i,j)=Qnet_ai(i,j)   ! net heat flux from ice to atmosphere
        END DO
      END DO
!
!  Open water case: set ice fluxes to zero.
!
      DO j = Jstr,Jend
        DO i = Istr,Iend
          IF (Si(i,j,linew,isAice).le.min_ai(ng)) THEN
            Fi(i,j,icIsst)=Fi(i,j,icT0mk)
            t2(i,j)=Fi(i,j,icT0mk)
            Si(i,j,linew,isTice)=-2.0_r8
#ifdef MASKING
            Fi(i,j,icIsst)=Fi(i,j,icIsst)*rmask(i,j)
            t2(i,j)=t2(i,j)*rmask(i,j)
            Si(i,j,linew,isTice)=Si(i,j,linew,isTice)*rmask(i,j)
# ifdef WET_DRY
            Fi(i,j,icIsst)=Fi(i,j,icIsst)*rmask_wet(i,j)
            t2(i,j)=t2(i,j)*rmask_wet(i,j)
            Si(i,j,linew,isTice)=Si(i,j,linew,isTice)*rmask_wet(i,j)
# endif
#endif
#ifdef ICESHELF
            IF (zice(i,j).ne.0.0_r8) THEN
              Fi(i,j,icIsst)=0.0_r8
              t2(i,j)=0.0_r8
              Si(i,j,linew,isTice)=0.0_r8
            END IF
#endif
            Qi2(i,j)=0.0_r8
            Qai(i,j)=0.0_r8
            Qio(i,j)=0.0_r8
            Si(i,j,linew,isHsno)=0.0_r8
            Si(i,j,linew,isHmel)=0.0_r8
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Suface water accumulation: ice melting.
!-----------------------------------------------------------------------
!
!  Set snow fall rate to value derived from precipitation rate.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ws(i,j)=MAX(snow(i,j), 0.0_r8)
        END DO
      END DO
!
!  Compute ice melt water thickness.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          tfrz=frln*sice(i,j)
          wsm(i,j)=0.0_r8
          Fi(i,j,icW_ai)=0.0_r8
          Fi(i,j,icW_ro)=0.0_r8
!
          IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
            cff=1.0_r8-brnfr(i,j)
            hfus1(i,j)=hfus*cff+                                        &
     &                 Fi(i,j,icIsst)*Cp_w-                             &
     &                 (cff*Cp_i+brnfr(i,j)*Cp_w)*Si(i,j,linew,isTice)
            Qai(i,j)=Qnet_ai(i,j)
            Qi2(i,j)=b2d(i,j)*(Si(i,j,linew,isTice)-Fi(i,j,icIsst))

            IF ((Si(i,j,linew,isHsno).le.eps).and.                      &
     &          (Si(i,j,linew,isHmel).le.eps)) THEN
              Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus1(i,j)*RhoSW)
            ELSE IF ((Si(i,j,linew,isHsno).le.eps).and.                 &
     &               (Si(i,j,linew,isHmel).gt.eps)) THEN
              Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus1(i,j)*1003.1_r8)
            ELSE
              Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus*SnowWetRho(ng))
            END IF

            IF ((Si(i,j,linew,isHsno).gt.eps).and.                      &
     &          (Fi(i,j,icIsst).ge.0.0_r8)) THEN
              Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)-                &
     &                             Si(i,j,linew,isAice)*                &
     &                             MAX(Qsur, 0.0_r8)*dtice(ng)
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             MAX(Qsur, 0.0_r8)*                   &
     &                             SnowWetRho(ng)/RhoSW*dtice(ng)
            ELSE IF ((Si(i,j,linew,isHmel).gt.eps).and.                 &
     &               (Fi(i,j,icIsst).le.tfrz)) THEN
              Fi(i,j,icW_ai)=MIN(Qsur, 0.0_r8)
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             MIN(Qsur, 0.0_r8)*dtice(ng)
            ELSE IF ((Si(i,j,linew,isHsno).le.eps).and.                 &
     &               (Si(i,j,linew,isHmel).ge.eps).and.                 &
     &               (Fi(i,j,icIsst).gt.tfrz)) THEN
              Fi(i,j,icW_ai)=MAX(Qsur, 0.0_r8)
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             MAX(Qsur, 0.0_r8)*dtice(ng)
            ELSE IF ((Si(i,j,linew,isHsno).lt.eps).and.                 &
     &               (Si(i,j,linew,isHmel).lt.eps).and.                 &
     &               (Fi(i,j,icIsst).gt.tfrz)) THEN
               Fi(i,j,icW_ai)=MAX(Qsur, 0.0_r8)
               Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+               &
     &                              Si(i,j,linew,isAice)*               &
     &                              MAX(Qsur, 0.0_r8)*dtice(ng)
            END IF

            IF (rain(i,j).le.0.0_r8) THEN
              Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             ws(i,j)*dtice(ng)
            ELSE IF ((Si(i,j,linew,isHsno).gt.0.0_r8).and.              &
     &               (Si(i,j,linew,isHmel).eq.0.0_r8)) THEN
              Si(i,j,linew,isHsno)=MAX(0.0_r8, Si(i,j,linew,isHsno)-    &
     &                             Si(i,j,linew,isAice)*rain(i,j)/      &
     &                             SnowDryRho(ng))
              Fi(i,j,icW_ai)=Fi(i,j,icW_ai)-                            &
     &                       2.0_r8*Si(i,j,linew,isAice)*               &
     &                       rain(i,j)/IceRho(ng)
            ELSE IF ((Si(i,j,linew,isHsno).gt.0.0_r8).and.              &
     &               (Si(i,j,linew,isHmel).gt.0.0_r8)) THEN
              Si(i,j,linew,isHsno)=MAX(0.0_r8, Si(i,j,linew,isHsno)-    &
     &                             0.5_r8*Si(i,j,linew,isAice)*         &
     &                             rain(i,j)/SnowDryRho(ng))
              Fi(i,j,icW_ai)=Fi(i,j,icW_ai)-                            &
     &                       0.5_r8*Si(i,j,linew,isAice)*               &
     &                       rain(i,j)/IceRho(ng)
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             0.5_r8*rain(i,j)/RhoSW*dtice(ng)
            ELSE IF (Si(i,j,linew,isHmel).gt.0.0_r8) THEN
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+                &
     &                             Si(i,j,linew,isAice)*                &
     &                             rain(i,j)/RhoSW*dtice(ng)
            ELSE
              Fi(i,j,icW_ai)=Fi(i,j,icW_ai)-                            &
     &                       Si(i,j,linew,isAice)*rain(i,j)/IceRho(ng)
            END IF
!
!  Limit the amount of surface water by the smaller of a max limit and
!  the ice thickness (SMD).
!
            IF (Si(i,j,linew,isHmel).gt.                                &
                MIN(max_hmelt(ng), Si(i,j,linew,isHice))) THEN
              Fi(i,j,icW_ro)=(Si(i,j,linew,isHmel)-                     &
     &                        MIN(max_hmelt(ng),                        &
     &                            Si(i,j,linew,isHice)))/dtice(ng)
              Si(i,j,linew,isHmel)=MIN(max_hmelt(ng),                   &
     &                                 Si(i,j,linew,isHice))
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Molecular sublayer under ice.
!-----------------------------------------------------------------------
!
!  Yaglom and Kader (1974) formulation for turbulent roughness length
!  scales "z0t" and "z0s".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          z0=MAX(z0ii*ice_thick(i,j), 0.01_r8)
          z0=MIN(z0, 0.1_r8)
          dztop=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          zdz0=dztop/z0
          zdz0=MAX(zdz0, 3.0_r8)
          rno=utau(i,j)*0.09_r8/nu
          termt=ykf*SQRT(rno)*prt**0.666667_r8
          terms=ykf*SQRT(rno)*prs**0.666667_r8
          cht(i,j)=utau(i,j)/(tpr*(LOG(zdz0)/kappa+termt))
          chs(i,j)=utau(i,j)/(tpr*(LOG(zdz0)/kappa+terms))
        END DO
      END DO
!
!  Temperature and salinity of molecular sublayer under ice.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          tfz=frln*salt_top(i,j)
          Fi(i,j,icW_ao)=0.0_r8
          Fi(i,j,icW_io)=0.0_r8
          xwai=MAX(0.0_r8, Fi(i,j,icW_ai))
          cff=1.0_r8-brnfr(i,j)
          hfus1(i,j)=hfus*cff+                                          &
     &               Fi(i,j,icT0mk)*Cp_w-                               &
     &               (cff*Cp_i+brnfr(i,j)*Cp_w)*Si(i,j,linew,isTice)
          IF (((temp_top(i,j).le.tfz).and.(Qnet_ao(i,j).gt.0.0_r8)).or. &
     &        ((temp_top(i,j).ge.tfz).and.(Qnet_ao(i,j).lt.0.0_r8).and. &
     &         (Si(i,j,linew,isAice).gt.0.0_r8))) THEN
            Fi(i,j,icW_ao)=Qnet_ao(i,j)/(hfus1(i,j)*RhoSW)
          END IF
          IF ((Si(i,j,linew,isAice).le.min_ai(ng)).or.                  &
     &        (Si(i,j,linew,isHice).le.min_hi(ng))) THEN
            Fi(i,j,icS0mk)=salt_top(i,j)
            Fi(i,j,icT0mk)=temp_top(i,j)
            Fi(i,j,icW_ai)=0.0_r8
            xtot=(1.0_r8-Si(i,j,linew,isAice))*Fi(i,j,icW_ao)
          ELSE                                          ! MK89 version
            Fi(i,j,icW_io)=(Qio(i,j)/RhoSW+                             &
     &                      Cp_w*cht(i,j)*(Fi(i,j,icT0mk)-              &
     &                                     temp_top(i,j)))/hfus1(i,j)
            xtot=Si(i,j,linew,isAice)*Fi(i,j,icW_io)+                   &
     &           (1.0_r8-Si(i,j,linew,isAice))*Fi(i,j,icW_ao)
!
!  Based on my reading of MK89, this calculation of "s0mk" does not
!  follow from the derivation. But it works quite well (SMD)!  Some
!  alternatives are commented below.
!
            Fi(i,j,icS0mk)=(chs(i,j)*salt_top(i,j)+                     &
     &                      (xwai-Fi(i,j,icW_io))*sice(i,j))/           &
     &                      (chs(i,j)+xwai+                             &
     &                       Fi(i,j,icW_ro)-Fi(i,j,icW_io))
!                                                               SMD s02
!           Fi(i,j,icS0mk)=(chs(i,j)*salt_top(i,j)+                     &
!    &                      (Si(i,j,linew,isAice)*Fi(i,j,icW_ro)-       &
!    &                       xtot)*sice(i,j))/                          &
!    &                      (chs(i,j)+
!    &                       Si(i,j,linew,isAice)*Fi(i,j,icW_ro)-xtot)
!
!                                                               SMD s03
!  Assume melt ponds are leaky, so replace "wro" with "xwai".
!
!           Fi(i,j,icS0mk)=(chs(i,j)*salt_top(i,j)+                     &
!    &                 (Si(i,j,linew,isAice)*xwai-xtot)*sice(i,j))/     &
!    &                (chs(i,j)+Si(i,j,linew,isAice)*xwai-xtot)
!
!                                                               SDM s04
!  Modify the original formulation by considering the balance
!  only over the ice covered portion of the grid cell such that
!  "wao" does not enter the expression.
!
!           Fi(i,j,icS0mk)=(chs(i,j)*salt_top(i,j)+                     &
!    &                      Si(i,j,linew,isAice)*                       &
!    &                      (xwai-Fi(i,j,icW_io))*sice(i,j))/           &
!    &                     (chs(i,j)+
!    &                      Si(i,j,linew,isAice)*xwai-Fi(i,j,icW_io))
!
!                                                               SMD s05
!  If we are to use what looks to be MK98 original formulation
!  it would have the P-E term in the denominator here as well.
!  Our PE term has the (1-ai) term factored into it already
!  this is used with our 'MKorig' formulation.
!
!           Fi(i,j,icS0mk)=(chs(i,j)*salt_top(i,j)+                     &
!    &                      (Si(i,j,linew,isAice)*xwai-                 &
!    &                       xtot)*sice(i,j))/                          &
!    &                     (chs(i,j)+Si(i,j,linew,isAice)*xwai-xtot+    &
!    &                      Si(i,j,linew,isAice)*stflx(i,j,isalt))
!
            Fi(i,j,icS0mk)=MAX(Fi(i,j,icS0mk), 0.0_r8)
            Fi(i,j,icS0mk)=MIN(Fi(i,j,icS0mk), 40.0_r8)
            Fi(i,j,icT0mk)=frln*Fi(i,j,icS0mk)
          END IF
!
!  Adjust surface heat flux.
!
          fac_shflx=1.0_r8
#ifdef ICESHELF
          IceCavity=zice(i,j).ne.0.0_r8
#else
          IceCavity=.FALSE.
#endif
          IF (.not.IceCavity) THEN
            IF(Si(i,j,linew,isAice).le.min_ai(ng)) THEN
               stflx(i,j,itemp)=Qnet_ao(i,j)*fac_shflx
            ELSE
#ifdef ICE_SHOREFAST
              hh=h(i,j)+Zt_avg1(i,j)
              clear=hh-0.9_r8*Si(i,j,liol,isHice)
              clear=MAX(clear, 0.0_r8)
              IF (clear.lt.1.5_r8) THEN
                fac_sf=MAX(clear-0.5_r8, 0.0_r8)/1.0_r8
              ELSE
                fac_sf=1.0_r8
              END IF
              stflx(i,j,itemp)=(1.0_r8-Si(i,j,linew,isAice))*           &
     &                         Qnet_ao(i,j)*fac_shflx+                  &
     &                         (Si(i,j,linew,isAice)*Qio(i,j)-          &
     &                          xtot*hfus1(i,j))*fac_sf
#else
              stflx(i,j,itemp)=(1.0_r8-Si(i,j,linew,isAice))*           &
     &                         Qnet_ao(i,j)+                            &
     &                         Si(i,j,linew,isAice)*Qio(i,j)-           &
     &                         xtot*hfus1(i,j)*RhoSW
#endif
            END IF
!
!  Change surface heat flux back to ROMS convention (Celsius m/s).
!
            stflx(i,j,itemp)=-stflx(i,j,itemp)*RhoCpr
#ifdef MASKING
            stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
#endif
!
!  Adjust surface freshwater flux.
!
#ifdef ICE_SHOREFAST
            cff=MIN(MAX(Fi(i,j,icS0mk), 0.0_r8), 60.0_r8)
            stflx(i,j,isalt)=stflx(i,j,isalt)-                          &
     &                       ((xtot-Si(i,j,linew,isAice)*xwai)*         &
     &                        (sice(i,j)-cff)+                          &
     &                        Si(i,j,linew,isAice)*                     &
     &                        Fi(i,j,icW_ro)*cff)*fac_sf
#else
            stflx(i,j,isalt)=stflx(i,j,isalt)+                          &
     &                       ((Si(i,j,linew,isAice)*                    &
     &                         (Fi(i,j,icW_io)-Fi(i,j,icW_ai))+         &
     &                         (1.0_r8-Si(i,j,linew,isAice))*           &
     &                         Fi(i,j,icW_ao)+                          &
     &                         Fi(i,j,icW_fr)))*                        &
     &                       (salt_top(i,j)-sice(i,j))-                 &
     &                       Si(i,j,linew,isAice)*                      &
     &                       (Fi(i,j,icW_ro)-xwai)*salt_top(i,j)
!
!  Fixed flux rate as a function of ice growth alone, MconsS case (SMD).
!
!           stflx(i,j,isalt)=stflx(i,j,isalt)+                          &
!    &                       (Si(i,j,linew,isAice)*                     &
!    &                        (Fi(i,j,icW_io)-Fi(i,j,icW_ai))+          &
!    &                        (1.0_r8-Si(i,j,linew,isAice))*            &
!    &                        Fi(i,j,icW_ao)+                           &
!    &                        Fi(i,j,icW_fr))*28.3_r8
!
!  If we want ice to have no effect on salinity (SMD):
!
!           IF ((Si(i,j,linew,isAice).gt.0.01_r8).and.                  &
!    &          (stflx(i,j,isalt).lt.0.0_r8)) THEN
!    &        stflx(i,j,isalt)=stflx(i,j,isalt)/                        &
!    &                         (1.0_r8-Si(i,j,linew,isAice))
!           END IF
!
!  Or alternatively we can include the precipitation over the ice (SMG).
!
!           stflx(i,j,isalt)=stflx(i,j,isalt)-                          &
!     &                      Si(i,j,linew,isAice)*                      &
!     &                      (Fi(i,j,icW_ro)-xwai)*salt_top(i,j)
#endif
#ifdef MASKING
            stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#endif
#ifdef WET_DRY
            stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#endif
!
!  Compute ice-ocean mass flux.
!
            Fi(i,j,icIOmf)=xtot-                                        &
     &                     Si(i,j,linew,isAice)*xwai-                   &
     &                     Si(i,j,linew,isAice)*Fi(i,j,icW_ro)+         &
     &                     Fi(i,j,icW_fr)
#ifdef MASKING
            Fi(i,j,icIOmf)=Fi(i,j,icIOmf)*rmask(i,j)
#endif
#ifdef WET_DRY
            Fi(i,j,icIOmf)=Fi(i,j,icIOmf)*rmask_wet(i,j)
#endif
          ELSE
            Fi(i,j,icIOmf)=0.0_r8
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Update ice properties.
!-----------------------------------------------------------------------
!
!  Track the amount of new ice produced thermodynamically to calculate
!  average ice age.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          mi_old=Si(i,j,linew,isHice)            ! old ice mass
          phi=3.0_r8
          IF (Fi(i,j,icW_ao).lt. 0.0_r8) phi=0.5_r8
          xmelt=MIN((Fi(i,j,icW_io)-Fi(i,j,icW_ai)), 0.0_r8)
          Si(i,j,linew,isHice)=Si(i,j,linew,isHice)+                    &
     &                         dtice(ng)*                               &
     &                         (Si(i,j,linew,isAice)*                   &
     &                          (Fi(i,j,icW_io)-Fi(i,j,icW_ai))+        &
     &                          (1.0_r8-Si(i,j,linew,isAice))*          &
     &                          Fi(i,j,icW_ao)+Fi(i,j,icW_fr))

          ai_tmp=Si(i,j,linew,isAice)            ! old ice concentration
          Si(i,j,linew,isAice)=Si(i,j,linew,isAice)+                    &
     &                         dtice(ng)*                               &
     &                         (1.0_r8-Si(i,j,linew,isAice))*           &
     &                         (phi*Fi(i,j,icW_ao)+Fi(i,j,icW_fr))
          Si(i,j,linew,isAice)= MIN(Si(i,j,linew,isAice), max_ai(ng))
          IF (Si(i,j,linew,isAice).lt.ai_tmp) THEN
            Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)*                  &
     &                           Si(i,j,linew,isAice)/MAX(ai_tmp, eps)
          END IF

#ifdef ICE_CONVSNOW
!
!  If snow base is below sea level, then raise the snow base to sea
!  level by converting some snow to ice (N.B. "hstar" is also weighted
!  by "isAice" like "isHsno" and "isHice").
!
          hstar=Si(i,j,linew,isHsno)-                                   &
     &          (RhoSW-IceRho(ng))*Si(i,j,linew,isHice)/SnowDryRho(ng)
          IF (hstar.gt.0.0_r8) THEN
            cff=hstar/RhoSW
            Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)-                  &
     &                           IceRho(ng)*cff
            Si(i,j,linew,isHice)=Si(i,j,linew,isHice)+                  &
     &                           SnowDryRho(ng)*cff
          END IF
#endif
#ifdef AICLM_NUDGING
          cff=AInudgcof(i,j)
          Si(i,j,linew,isAice)=Si(i,j,linew,isAice)+                    &
     &                         dtice(ng)*cff*                           &
     &                         (aiclm(i,j)-Si(i,j,linew,isAice))
          Si(i,j,linew,isHice)=Si(i,j,linew,isHice)+                    &
     &                         dtice(ng)*cff*                           &
     &                         (hiclm(i,j)-Si(i,j,linew,isHice))
#endif
!
!  Determine age of the sea ice. Any new ice production reduces the
!  overall age of the ice parcel.
!
          IF ((Si(i,j,linew,isIage).le.0.0_r8).and.                     &
     &        (Si(i,j,linew,isHice).gt.min_hi(ng))) THEN        ! new
            Si(i,j,linew,isIage)=dtice(ng)*sec2day
          ELSE IF((Si(i,j,linew,isIage).gt.0.0_r8).and.                 &
     &            (Si(i,j,linew,isHice).gt.min_hi(ng))) THEN    ! older
            delta_mi=MIN(MAX(Si(i,j,linew,isHice)-mi_old, 0.0_r8)/      &
     &                       Si(i,j,linew,isHice), 1.0_r8)
            Si(i,j,linew,isIage)=Si(i,j,linew,isIage)+                  &
     &                           dtice(ng)*sec2day-                     &
     &                           Si(i,j,linew,isIage)*delta_mi
          ELSE                                                  ! melted
            Si(i,j,linew,isIage)=0.0_r8
          ENDIF

#ifdef MASKING
          Si(i,j,linew,isAice)=Si(i,j,linew,isAice)*rmask(i,j)
          Si(i,j,linew,isHice)=Si(i,j,linew,isHice)*rmask(i,j)
#endif
#ifdef WET_DRY
!         Si(i,j,linew,isAice)=Si(i,j,linew,isAice)*rmask_wet(i,j)
!         Si(i,j,linew,isHice)=Si(i,j,linew,isHice)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
            Si(i,j,linew,isAice)=0.0_r8
            Si(i,j,linew,isHice)=0.0_r8
          END IF
#endif
        END DO
      END DO
!
!  Limit the values.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Si(i,j,linew,isAice)=MIN(Si(i,j,linew,isAice), max_ai(ng))
          Si(i,j,linew,isAice)=MAX(Si(i,j,linew,isAice), 0.0_r8)
          Si(i,j,linew,isHice)=MAX(Si(i,j,linew,isHice), 0.0_r8)
          Si(i,j,linew,isHsno)=MAX(Si(i,j,linew,isHsno), 0.0_r8)
          Si(i,j,linew,isHmel)=MAX(Si(i,j,linew,isHmel), 0.0_r8)
          Si(i,j,linew,isTice)=MAX(Si(i,j,linew,isTice), -70.0_r8)
          IF (Si(i,j,linew,isHice).le.0.0_r8)                           &
     &      Si(i,j,linew,isAice)=0.0_r8
          IF (Si(i,j,linew,isAice).le.0.0_r8)                           &
     &      Si(i,j,linew,isHice)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icIsst))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icQcon))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icQrhs))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))

      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsAice),                             &
     &                    LBC(:,ibICE(isAice),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsHice),                             &
     &                    LBC(:,ibICE(isHice),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsHsno),                             &
     &                    LBC(:,ibICE(isHsno),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsHmel),                             &
     &                    LBC(:,ibICE(isHmel),ng))

      CALL ice_bc2d_tile (ng, tile, model, isIage,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsIage),                             &
     &                    LBC(:,ibICE(isIage),ng))

      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,IsUice),                             &
     &                    Si(:,:,:,IsVice),                             &
     &                    Si(:,:,:,IsHice),                             &
     &                    Si(:,:,:,IsTice),                             &
     &                    Si(:,:,:,IsEnth))
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHage))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHmel))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHsno))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isAice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isIage))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isEnth))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isTice))
      END IF

#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isHage),                         &
     &                    Si(:,:,linew,isHice),                         &
     &                    Si(:,:,linew,isHmel),                         &
     &                    Si(:,:,linew,isHsno))

      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isAice),                         &
     &                    Si(:,:,linew,isIage),                         &
     &                    Si(:,:,linew,isEnth),                         &
     &                    Si(:,:,linew,isTice))
#endif
!
      RETURN
      END SUBROUTINE ice_thermo_tile
!
      END MODULE ice_thermo_mod
