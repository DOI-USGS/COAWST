      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Nemuro Lower Trophic Level Ecosystem Model.                         !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields. Currently, the ecosystem has     !
!  the following functional compartments:                              !
!                                                                      !
!    iSphy     small phytoplankton biomass, nanophytoplankton          !
!    iLphy     large phytoplankton biomass, diatoms                    !
!    iSzoo     small zooplankton biomass, microzooplankton (ciliates)  !
!    iLzoo     large zooplankton biomass, mesozooplankton (copepods)   !
!    iPzoo     predator zooplankton biomass (euphausiids, etc)         !
!    iNO3_     nitrate concentration, NO3                              !
!    iNH4_     ammonium concentration, NH4                             !
!    iPON_     particulate organic nitrogen                            !
!    iDON_     dissolved organic nitrogen                              !
!    iSiOH     silicate concentration, Si(OH)4 (silicic acid)          !
!    iopal     particulate organic silica                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Kishi, M. J., et  all, 2007: Nemuro - a lower trophic level       !
!      model for the North Pacific marine ecosystem,  Ecological       !
!      Modelling, 202, 12-25.                                          !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
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
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
#ifdef IRON_LIMIT
     &                   GRID(ng) % h,                                  &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#ifdef NEMURO_SED1
     &                   OCEAN(ng) % PONsed,                            &
     &                   OCEAN(ng) % OPALsed,                           &
     &                   OCEAN(ng) % DENITsed,                          &
     &                   OCEAN(ng) % PON_burial,                        &
     &                   OCEAN(ng) % OPAL_burial,                       &
#endif
#ifdef PRIMARY_PROD
     &                   OCEAN(ng) % Bio_NPP,                           &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
#ifdef IRON_LIMIT
     &                         h,                                       &
#endif
     &                         Hz, z_r, z_w,                            &
     &                         srflx,                                   &
#ifdef NEMURO_SED1
     &                         PONsed, OPALsed, DENITsed,               &
     &                         PON_burial, OPAL_burial,                 &
#endif
#ifdef PRIMARY_PROD
     &                         Bio_NPP,                                 &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
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
# endif
#ifdef IRON_LIMIT
      real(r8), intent(in) :: h(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# ifdef NEMURO_SED1
      real(r8), intent(inout) :: PONsed(LBi:,LBj:)
      real(r8), intent(inout) :: OPALsed(LBi:,LBj:)
      real(r8), intent(inout) :: DENITsed(LBi:,LBj:)
      real(r8), intent(inout) :: PON_burial(LBi:,LBj:)
      real(r8), intent(inout) :: OPAL_burial(LBi:,LBj:)
# endif
# ifdef PRIMARY_PROD
      real(r8), intent(out) :: Bio_NPP(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef IRON_LIMIT
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# ifdef NEMURO_SED1
      real(r8), intent(inout) :: PONsed(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: OPALsed(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: DENITsed(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: PON_burial(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: OPAL_burial(LBi:UBi,LBj:UBj)
# endif
# ifdef PRIMARY_PROD
      real(r8), intent(out) :: Bio_NPP(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 2

      integer :: Iter, ibio, indx, isink, itime, itrc, iTrcMax
      integer :: i, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: MinVal = 1.0e-6_r8

      real(r8) :: AttL, AttS, IrrL, IrrS, KappaL, KappaS
      real(r8) :: dtdays, dz
      real(r8) :: GppAPS, GppAPL, GppNPS, GppNPL, GppPS, GppPL
      real(r8) :: GraPL2ZS, GraPL2ZL, GraPL2ZP, GraPS2ZL, GraPS2ZS
      real(r8) :: GraZL2ZP, GraZS2ZL, GraZS2ZP
!      real(r8) :: GraZS2PS, GraZL2PS, GraZL2PL, GraZL2ZS
!      real(r8) :: GraZP2PL, GraZP2ZS, GraZP2ZL
      real(r8) :: EgeZL, EgeZP, EgeZS
      real(r8) :: ExcPL, ExcPS, ExcZL, ExcZP, ExcZS
      real(r8) :: MorPL, MorPS
      real(r8) :: ResPL, ResPS
      real(r8) :: RnewL, RnewS
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8
      real(r8) :: fac, fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: Flimit
#ifdef IRON_LIMIT
      real(r8) :: FNratio, FCratio, FCratioE
      real(r8) :: cffFe, Fndgcf
      real(r8) :: h_max, Fe_min, Fe_max, Fe_rel, SiN_min, SiN_max
#endif

      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur

      real(r8), dimension(NT(ng),2) :: BioTrc
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: LightL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: LightS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      real(r8), dimension(IminS:ImaxS,N(ng)) :: LPSiN
#ifdef PRIMARY_PROD
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NPP_slice
#endif
      real(r8), dimension(IminS:ImaxS) :: frac_buried

#include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping size (days) according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPON_                 ! particulate organic nitrogen
      idsink(2)=iopal                 ! particulate organic silica
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=setVPON(ng)             ! particulate organic nitrogen
      Wbio(2)=setVOpal(ng)            ! particulate organic silica
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
#ifdef PRIMARY_PROD
        DO k=1,N(ng)
          DO i=Istr,Iend
            NPP_slice(i,k)=0.0_r8
          END DO
        END DO
#endif
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
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
!#ifdef ANA_BIOSWRAD
!          cff=2.0_r8*3.14159_r8*time(ng)/86400.0_r8/360.0_r8
!          PARsur(i)=PARfrac(ng)*200.0_r8*0.5_r8*(1.0_r8-cos(cff))
#ifdef CONST_PAR
          PARsur(i)=PARfrac(ng)*200.0_r8
#else
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#endif
        END DO
!
#if defined IRON_LIMIT && defined IRON_RELAX
!  Relaxation of dissolved iron to climatology
        DO k=1,N(ng)
          DO i=Istr,Iend
!  Set concentration and depth parameters for FeD climatology
!  Fe concentration in (micromol-Fe/m3, or nM-Fe)
            h_max = 200.0_r8
            Fe_max = 2.0_r8
!  Set nudging time scales to 5 days
            Fe_rel = 5.0_r8
            Fndgcf = 1.0_r8/(Fe_rel*86400.0_r8)
!  Relaxation for depths < h_max to simulate Fe input at coast
            IF (h(i,j).le.h_max) THEN
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       dt(ng)*Fndgcf*(Fe_max-Bio(i,k,iFeD_))
            ELSE IF (h(i,j).gt.h_max .and. Bio(i,k,iFeD_).lt.0.1_r8) THEN
                    Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       dt(ng)*Fndgcf*(0.1_r8-Bio(i,k,iFeD_))
            END IF
          END DO
        END DO
#endif
!
#if defined IRON_LIMIT && defined IRON_RSIN
!  Variable Si:N ratio for dimatoms based on dissolved iron concentration
!  Si:N varies between 1:1 (high iron) and 3:1 (low iron)
        DO k=1,N(ng)
          DO i=Istr,Iend
            Fe_min = 0.05_r8
            Fe_max = 2.0_r8
            SiN_min = 1.0_r8
            SiN_max = 3.0_r8
            cffFe=max(0.0_r8,min(1.0_r8,(Bio(i,k,iFeD_)-Fe_min)/        &
     &              (Fe_max-Fe_min)))
            LPSiN(i,k)=SiN_max-cffFe*(SiN_max-SiN_min)
          END DO
        END DO
#else
        DO k=1,N(ng)
          DO i=Istr,Iend
            LPSiN(i,k)=RSiN(ng)
          END DO
        END DO
#endif
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
!  in a fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!  In the implicit algorithm, we have for example (N: nutrient,
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
!  Compute light attenuation as function of depth.
!
          cff1=1.0/VmaxS(ng)
          cff2=1.0/VmaxL(ng)
          DO i=Istr,Iend
            AttS=PARsur(i)
            AttL=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell using the
!  Platt et al. (1980) photoinhibition formulation. Here, AttSW is
!  the light attenuation due to seawater and AttPS and AttPL is the
!  attenuation due to Small and Large Phytoplankton (self-shading
!  coefficient).
!
                dz=0.5_r8*(z_w(i,j,k)-z_w(i,j,k-1))
                kappaS=AttSW(ng)+                                       &
     &                 AttPS(ng)*(Bio(i,k,iSphy)+Bio(i,k,iLphy))
                kappaL=AttSW(ng)+                                       &
     &                 AttPL(ng)*(Bio(i,k,iSphy)+Bio(i,k,iLphy))
                IrrS=EXP(-kappaS*dz)
                IrrL=EXP(-kappaL*dz)
                AttS=AttS*IrrS
                AttL=AttL*IrrL
                LightS(i,k)=(1.0_r8-EXP(-alphaPS(ng)*AttS*cff1))*       &
     &                      EXP(-betaPS(ng)*AttS*cff1)
                LightL(i,k)=(1.0_r8-EXP(-alphaPL(ng)*AttL*cff2))*       &
     &                      EXP(-betaPL(ng)*AttL*cff2)
!
!  Attenuate the light to the bottom of the grid cell.
!
                AttS=AttS*IrrS
                AttL=AttL*IrrL
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                LightS(i,k)=0.0_r8
                LightL(i,k)=0.0_r8
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton primary productivity.
!-----------------------------------------------------------------------
!
!  Gross primary production of Small Phytoplankton consisting of
!  nutrient uptake (NO3 and NH4) terms, temperature-dependend term,
!  and light limitation term. The Michaelis-Menten curve is used to
!  describe the change in uptake rate as a function of nutrient
!  concentration.
!
          cff=dtdays*VmaxS(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
#ifdef IRON_LIMIT
! Small phytoplankton growth reduction factor due to iron limitation
!
! Current Fe:N ratio [umol-Fe/mmol-N]
              FNratio=Bio(i,k,iFeSp)/MAX(MinVal,Bio(i,k,iSphy))
! Current F:C ratio [umol-Fe/mol-C]
! (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1e3 mmol-C/mol-C),
              FCratio=FNratio*(16.0_r8/106.0_r8)*1.0e3_r8
! Empirical FCratio
              FCratioE= B_Fe(ng)*Bio(i,k,iFeD_)**A_Fe(ng)
! Phytoplankton growth reduction factor due to iron limitation
! based on F:C ratio
              Flimit = FCratio**2.0_r8/                                 &
     &                 (FCratio**2.0_r8+SK_FeC(ng)**2.0_r8)
!JF              Flimit = FCratio/(FCratio+SK_FeC(ng))
#else
              Flimit = 1.0_r8
#endif
!
!  Small Phytoplankton gross primary productivity, GppPS.
!
              cff1=cff*EXP(KGppS(ng)*Bio(i,k,itemp))*LightS(i,k)*       &
     &             Bio(i,k,iSphy)
              cff2=EXP(-PusaIS(ng)*Bio(i,k,iNH4_))/                     &
     &             (KNO3S(ng)+Bio(i,k,iNO3_))
              cff3=1.0_r8/(KNH4S(ng)+Bio(i,k,iNH4_))
              cff4=cff2*Bio(i,k,iNO3_)
              cff5=cff3*Bio(i,k,iNH4_)
              cff6=Flimit/MAX(MinVal,cff4+cff5)
              cff4=cff1*cff2*MIN(1.0_r8,cff6)        ! Fe lim on N03
              cff5=cff1*cff3*MIN(1.0_r8,cff6)        ! Fe lim on NH4
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
              GppNPS=Bio(i,k,iNO3_)*cff4
              GppAPS=Bio(i,k,iNH4_)*cff5
              GppPS=GppNPS+GppAPS
              Bio(i,k,iSphy)=Bio(i,k,iSphy)+GppPS
#ifdef PRIMARY_PROD
              NPP_slice(i,k)=NPP_slice(i,k)+GppPS
#endif
!
#ifdef IRON_LIMIT
! Iron uptake proportional to growth
              cffFe=GppPS*FNratio/MAX(MinVal,Bio(i,k,iFeD_))
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)+                            &
     &                       Bio(i,k,iFeD_)*cffFe
!
! Iron uptake to reach appropriate Fe:C ratio
              cffFe=dtdays*(FCratioE-FCratio)/T_fe(ng)
              cffFe=Bio(i,k,iSphy)*cffFe*(106.0_r8/16.0_r8)*1.0e-3_r8
              IF (cffFe.ge.0.0_r8) THEN
                cffFe=cffFe/MAX(MinVal,Bio(i,k,iFeD_))
                Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
                Bio(i,k,iFeSp)=Bio(i,k,iFeSp)+                          &
     &                         Bio(i,k,iFeD_)*cffFe
              ELSE
                cffFe=-cffFe/MAX(MinVal,Bio(i,k,iFeSp))
                Bio(i,k,iFeSp)=Bio(i,k,iFeSp)/(1.0_r8+cffFe)
                Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                          &
     &                         Bio(i,k,iFeSp)*cffFe
              END IF
#endif
!
!  Small Phytoplankton respiration rate, ResPS, assumed to be
!  proportional to biomass. Use ratio of NO3 uptake to total update
!  (NO3 + NH4) to compute respiration contributions.
!
              RnewS=GppNPS/MAX(MinVal,GppPS)
              cff4=dtdays*ResPS0(ng)*EXP(KResPS(ng)*Bio(i,k,itemp))
              Bio(i,k,iSphy)=Bio(i,k,iSphy)/(1.0_r8+cff4)
              ResPS=Bio(i,k,iSphy)*cff4
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+ResPS*RnewS
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+ResPS*(1.0_r8-RnewS)
#ifdef PRIMARY_PROD
              NPP_slice(i,k)=NPP_slice(i,k)-ResPS
#endif
!
#ifdef IRON_LIMIT
! Small phytoplankton respiration Fe terms
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)/(1.0_r8+cff4)
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       Bio(i,k,iFeSp)*cff4*FeRR(ng)
#endif
!
!  Small Phytoplankton extracellular excrection rate, ExcPS, assumed to
!  be proportional to production.
!
              ExcPS=GppPS*GammaS(ng)
              Bio(i,k,iSphy)=Bio(i,k,iSphy)-ExcPS
              Bio(i,k,iDON_)=Bio(i,k,iDON_)+ExcPS
!
#ifdef IRON_LIMIT
! Small phytoplankton excretion Fe terms
              cffFe=ExcPS*Bio(i,k,iFeSp)/MAX(MinVal,Bio(i,k,iSphy))
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
            END DO
          END DO
!
!  Gross primary production of Large Phytoplankton consisting of
!  nutrient uptake (NO3, NH4, Silicate) terms, temperature-dependend
!  term, and light limitation term. Notice that there is a silicate
!  limiting term (which complicates the implicit algorithm). Again,
!  the Michaelis-Menten curve is used to describe the change in
!  uptake rate as a function of nutrient concentration.
!
          cff=dtdays*VmaxL(ng)
!          fac1=1.0/RSiN(ng)
          fac2=dtdays*ResPL0(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
#ifdef IRON_LIMIT
! Large phytoplankton growth reduction factor due to iron limitation
!
! Current Fe:N ratio [umol-Fe/mmol-N]
              FNratio=Bio(i,k,iFeLp)/MAX(MinVal,Bio(i,k,iLphy))
! Current F:C ratio [umol-Fe/mol-C]
! (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1e3 mmol-C/mol-C),
              FCratio=FNratio*(16.0_r8/106.0_r8)*1.0e3_r8
! Empirical FCratio
              FCratioE= B_Fe(ng)*Bio(i,k,iFeD_)**A_Fe(ng)
! Phytoplankton growth reduction factor due to Iron limitation
! based on F:C ratio
              Flimit = FCratio**2.0_r8/                                 &
     &                 (FCratio**2.0_r8+LK_FeC(ng)**2.0_r8)
!JF              Flimit = FCratio/(FCratio+LK_FeC(ng))
#else
              Flimit = 1.0_r8
#endif
!
!  Large Phytoplankton gross primary productivity, GppPL.
!
              cff1=cff*EXP(KGppL(ng)*Bio(i,k,itemp))*LightL(i,k)*       &
     &             Bio(i,k,iLphy)
              cff2=EXP(-PusaIL(ng)*Bio(i,k,iNH4_))/                     &
     &             (KNO3L(ng)+Bio(i,k,iNO3_))
              cff3=1.0_r8/(KNH4L(ng)+Bio(i,k,iNH4_))
              cff4=cff2*Bio(i,k,iNO3_)
              cff5=cff3*Bio(i,k,iNH4_)
              cff6=Bio(i,k,iSiOH)/(KSiL(ng)+Bio(i,k,iSiOH))
              cff7=cff6/MAX(MinVal,cff4+cff5)
              cff8=Flimit/MAX(MinVal,cff4+cff5)
              cff4=cff1*cff2*MIN(1.0_r8,MIN(cff7,cff8))   ! Si/Fe lim on N03
              cff5=cff1*cff3*MIN(1.0_r8,MIN(cff7,cff8))   ! Si/Fe lim on NH4
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
              GppNPL=Bio(i,k,iNO3_)*cff4
              GppAPL=Bio(i,k,iNH4_)*cff5
              GppPL=GppNPL+GppAPL
              Bio(i,k,iLphy)=Bio(i,k,iLphy)+GppPL
              Bio(i,k,iSiOH)=Bio(i,k,iSiOH)-GppPL*LPSiN(i,k)
#ifdef PRIMARY_PROD
              NPP_slice(i,k)=NPP_slice(i,k)+GppPL
#endif
!
#ifdef IRON_LIMIT
! Iron Uptake proportional to growth
              cffFe=GppPL*FNratio/MAX(MinVal,Bio(i,k,iFeD_))
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)+                            &
     &                       Bio(i,k,iFeD_)*cffFe
!
! Iron uptake to reach appropriate Fe:C ratio
              cffFe=dtdays*(FCratioE-FCratio)/T_fe(ng)
              cffFe=Bio(i,k,iLphy)*cffFe*(106.0_r8/16.0_r8)*1.0e-3_r8
              IF (cffFe.ge.0.0_r8) THEN
                cffFe=cffFe/MAX(MinVal,Bio(i,k,iFeD_))
                Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
                Bio(i,k,iFeLp)=Bio(i,k,iFeLp)+                          &
     &                         Bio(i,k,iFeD_)*cffFe
              ELSE
                cffFe=-cffFe/MAX(MinVal,Bio(i,k,iFeLp))
                Bio(i,k,iFeLp)=Bio(i,k,iFeLp)/(1.0_r8+cffFe)
                Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                          &
     &                         Bio(i,k,iFeLp)*cffFe
              END IF
#endif
!
!  Large Phytoplankton respiration rate, ResPL, assumed to be
!  proportional to biomass. Use ratio of NO3 uptake to total update
!  (NO3 + NH4) to compute respiration contributions. Use Si:N ratio to
!  compute SiOH4 contribution.
!
              RnewL=GppNPL/MAX(MinVal,GppPL)
              cff7=fac2*EXP(KResPL(ng)*Bio(i,k,itemp))
              Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff7)
              ResPL=Bio(i,k,iLphy)*cff7
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+ResPL*RnewL
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+ResPL*(1.0_r8-RnewL)
              Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+ResPL*LPSiN(i,k)
#ifdef PRIMARY_PROD
              NPP_slice(i,k)=NPP_slice(i,k)-ResPL
#endif
!
#ifdef IRON_LIMIT
! Large Phytoplankton respiration Fe terms
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)/(1.0_r8+cff5)
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       Bio(i,k,iFeLp)*cff5*FeRR(ng)
#endif
!
!  Large Phytoplankton extracellular excrection rate, ExcPL, assumed to
!  be proportional to production.
!
              ExcPL=GppPL*GammaL(ng)
              Bio(i,k,iLphy)=Bio(i,k,iLphy)-ExcPL
              Bio(i,k,iDON_)=Bio(i,k,iDON_)+ExcPL
              Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+ExcPL*LPSiN(i,k)
!
#ifdef IRON_LIMIT
! Large Phytoplankton excretion Fe terms
              cffFe=ExcPL*Bio(i,k,iFeLp)/MAX(MinVal,Bio(i,k,iLphy))
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton mortality to particulate organic nitrogen.
!-----------------------------------------------------------------------
!
          fac1=dtdays*MorPS0(ng)
          fac2=dtdays*MorPL0(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*Bio(i,k,iSphy)*EXP(KMorPS(ng)*Bio(i,k,itemp))
              cff2=fac2*Bio(i,k,iLphy)*EXP(KMorPL(ng)*Bio(i,k,itemp))
              Bio(i,k,iSphy)=Bio(i,k,iSphy)/(1.0_r8+cff1)
              Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff2)
              MorPS=Bio(i,k,iSphy)*cff1
              MorPL=Bio(i,k,iLphy)*cff2
              Bio(i,k,iPON_)=Bio(i,k,iPON_)+MorPS+MorPL
              Bio(i,k,iopal)=Bio(i,k,iopal)+MorPL*LPSiN(i,k)

#ifdef IRON_LIMIT
! Phytoplankton Mortality Fe terms
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)/(1.0_r8+cff1)
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       Bio(i,k,iFeSp)*cff1*FeRR(ng)
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)/(1.0_r8+cff2)
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       Bio(i,k,iFeLp)*cff2*FeRR(ng)
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton grazing, egestion and excretion.
!-----------------------------------------------------------------------

#if defined IVLEV_EXPLICIT
!
!  The rate of grazing by the zooplankton is modeled using an Ivlev
!  equation with a feeding threshold using an explicit (non-conserving)
!  algorithm to the original formulation. Notice that term is forced to
!  be positive using a MAX function.
!
#elif defined HOLLING_GRAZING
!
!  The rate of grazing by the zooplankton is modeled using a Holling-
!  type s-shaped curve. It is known to be numerically more stable and
!  allows an implicit discretization.
!
!    P(new) = P(old) - dt*mu*[P(old)/(Kp + P(old)^2)]*P(new)*Z(old)
!
!  The implicit grazing term is then:
!
!    P(new) = P(old) / (1 + G)
!
!  were the grazing rate, G, is:
!
!    G = dt * mu * [P(old)/(Kp + P(old)^2)] * Z(old)
!
#else
# define IVLEV_IMPLICIT
!
!  The rate of grazing by the zooplankton is modeled using an Ivlev
!  equation with a feeding threshold. An implicit algorithm is
!  achieved by multiplying grazing term by a unity factor, alpha:
!
!    alpha = 1 = P(new)/(P(old)+deltaP)
!
!  where
!
!    deltaP = P(new) - P(old) = - dt*mu*[1-EXP(lambda*P(old))]
!
!  The factor alpha can be approximated using Taylor series to:
!
!    alpha = [P(new) / P(old)] * [1 - deltaP/P(old)]
!
!  The discretized grazing term is then:
!
!    P(new) = P(old) - dt*mu*[1-EXP(lambda*P)] * alpha * Z(old)
!
!  Which can be approximated with an implicit algorithm to:
!
!    P(new) = P(old) / (1 + G)
!
!  were the grazing rate, G, is:
!
!    G = [1 + P(old)/(dt*mu*(1 - EXP(lambda*P(old))))] * Z(old)
!
#endif
          fac1=dtdays*GRmaxSps(ng)
          fac2=dtdays*GRmaxSpl(ng)
          fac3=dtdays*GRmaxLps(ng)
          fac4=dtdays*GRmaxLpl(ng)
          fac5=dtdays*GRmaxLzs(ng)
          fac6=dtdays*GRmaxPpl(ng)
          fac7=dtdays*GRmaxPzs(ng)
          fac8=dtdays*GRmaxPzl(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  Temperature-dependent term (Q10).
!
              cff1=EXP(KGraS(ng)*Bio(i,k,itemp))
              cff2=EXP(KGraL(ng)*Bio(i,k,itemp))
              cff3=EXP(KGraP(ng)*Bio(i,k,itemp))
!
!  Small Zooplankton grazing on Small Phytoplankton, GraPS2ZS.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamS(ng)*(PS2ZSstar(ng)-Bio(i,k,iSphy)))
              GraPS2ZS=fac1*cff1*MAX(0.0_r8,cff4)*Bio(i,k,iSzoo)
              Bio(i,k,iSphy)=Bio(i,k,iSphy)-GraPS2ZS
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)+GraPS2ZS
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KPS2ZS(ng)+Bio(i,k,iSphy)*Bio(i,k,iSphy))
              cff=fac1*cff1*cff4*Bio(i,k,iSzoo)*Bio(i,k,iSphy)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamS(ng)*(PS2ZSstar(ng)-Bio(i,k,iSphy)))
              cff5=1.0_r8/(fac1*cff4)
              cff=(1.0_r8+Bio(i,k,iSphy)*cff5)*cff1*Bio(i,k,iSzoo)
# endif
              Bio(i,k,iSphy)=Bio(i,k,iSphy)/(1.0_r8+cff)
              GraPS2ZS=cff*Bio(i,k,iSphy)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)+GraPS2ZS
#endif
#ifdef IRON_LIMIT
!     Grazing Fe terms
              cffFe=GraPS2ZS*Bio(i,k,iFeSp)/MAX(MinVal,Bio(i,k,iSphy))
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
!
!  Small Zooplankton grazing on Large Phytoplankton, GraPL2ZS.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PL2ZSstar(ng)-Bio(i,k,iLphy)))
              GraPL2ZS=fac2*cff2*MAX(0.0_r8,cff4)*Bio(i,k,iSzoo)
              Bio(i,k,iLphy)=Bio(i,k,iLphy)-GraPL2ZS
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)+GraPL2ZS
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KPL2ZS(ng)+Bio(i,k,iLphy)*Bio(i,k,iLphy))
              cff=fac2*cff2*cff4*Bio(i,k,iSzoo)*Bio(i,k,iLphy)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PL2ZSstar(ng)-Bio(i,k,iLphy)))
              cff5=1.0_r8/(fac2*cff4)
              cff=(1.0_r8+Bio(i,k,iLphy)*cff5)*cff2*Bio(i,k,iSzoo)
# endif
              Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff)
              GraPL2ZS=cff*Bio(i,k,iLphy)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)+GraPL2ZS
#endif
#ifdef IRON_LIMIT
!     Grazing Fe terms
              cffFe=GraPL2ZS*Bio(i,k,iFeLp)/MAX(MinVal,Bio(i,k,iLphy))
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
!
!  Large Zooplankton grazing on Small Phytoplankton, GraPS2ZL.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PS2ZLstar(ng)-Bio(i,k,iSphy)))
              GraPS2ZL=fac3*cff2*MAX(0.0_r8,cff4)*Bio(i,k,iLzoo)
              Bio(i,k,iSphy)=Bio(i,k,iSphy)-GraPS2ZL
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraPS2ZL
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KPS2ZL(ng)+Bio(i,k,iSphy)*Bio(i,k,iSphy))
              cff=fac3*cff2*cff4*Bio(i,k,iLzoo)*Bio(i,k,iSphy)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PS2ZLstar(ng)-Bio(i,k,iSphy)))
              cff5=1.0_r8/(fac3*cff4)
              cff=(1.0_r8+Bio(i,k,iSphy)*cff5)*cff2*Bio(i,k,iLzoo)
# endif
              Bio(i,k,iSphy)=Bio(i,k,iSphy)/(1.0_r8+cff)
              GraPS2ZL=cff*Bio(i,k,iSphy)
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraPS2ZL
#endif
#ifdef IRON_LIMIT
!     Grazing Fe terms
              cffFe=GraPS2ZL*Bio(i,k,iFeSp)/MAX(MinVal,Bio(i,k,iSphy))
              Bio(i,k,iFeSp)=Bio(i,k,iFeSp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
!
!  Large Zooplankton grazing on Large Phytoplankton, GraPL2ZL.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PL2ZLstar(ng)-Bio(i,k,iLphy)))
              GraPL2ZL=fac4*cff2*MAX(0.0_r8,cff4)*Bio(i,k,iLzoo)
              Bio(i,k,iLphy)=Bio(i,k,iLphy)-GraPL2ZL
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraPL2ZL
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KPL2ZL(ng)+Bio(i,k,iLphy)*Bio(i,k,iLphy))
              cff=fac4*cff2*cff4*Bio(i,k,iLzoo)*Bio(i,k,iLphy)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(PL2ZLstar(ng)-Bio(i,k,iLphy)))
              cff5=1.0_r8/(fac4*cff4)
              cff=(1.0_r8+Bio(i,k,iLphy)*cff5)*cff2*Bio(i,k,iLzoo)
# endif
              Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff)
              GraPL2ZL=cff*Bio(i,k,iLphy)
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraPL2ZL
#endif
#ifdef IRON_LIMIT
!     Grazing Fe terms
              cffFe=GraPL2ZL*Bio(i,k,iFeLp)/MAX(MinVal,Bio(i,k,iLphy))
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
!
!  Large Zooplankton grazing on Small Zooplankton, GraZS2ZL.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(ZS2ZLstar(ng)-Bio(i,k,iSzoo)))
              GraZS2ZL=fac5*cff2*MAX(0.0_r8,cff4)*Bio(i,k,iLzoo)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)-GraZS2ZL
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraZS2ZL
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KZS2ZL(ng)+Bio(i,k,iSzoo)*Bio(i,k,iSzoo))
              cff=fac5*cff2*cff4*Bio(i,k,iLzoo)*Bio(i,k,iSzoo)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamL(ng)*(ZS2ZLstar(ng)-Bio(i,k,iSzoo)))
              cff5=1.0_r8/(fac5*cff4)
              cff=(1.0_r8+Bio(i,k,iSphy)*cff5)*cff2*Bio(i,k,iLzoo)
# endif
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)/(1.0_r8+cff)
              GraZS2ZL=cff*Bio(i,k,iSzoo)
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+GraZS2ZL
#endif
!
!  Predactor Zooplankton grazing on Large Phytoplankton, GraPL2ZP.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(PL2ZPstar(ng)-Bio(i,k,iLphy)))
              cff5=EXP(-PusaiPL(ng)*(Bio(i,k,iLzoo)+Bio(i,k,iSzoo)))
              GraPL2ZP=fac6*cff3*cff5*MAX(0.0_r8,cff4)*Bio(i,k,iPzoo)
              Bio(i,k,iLphy)=Bio(i,k,iLphy)-GraPL2ZP
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraPL2ZP
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KPL2ZP(ng)+Bio(i,k,iLphy)*Bio(i,k,iLphy))
              cff5=EXP(-PusaiPL(ng)*(Bio(i,k,iLzoo)+Bio(i,k,iSzoo)))
              cff=fac6*cff3*cff4*cff5*Bio(i,k,iPzoo)*Bio(i,k,iLphy)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(PL2ZPstar(ng)-Bio(i,k,iLphy)))
              cff5=EXP(-PusaiPL(ng)*(Bio(i,k,iLzoo)+Bio(i,k,iSzoo)))
              cff6=1.0_r8/(fac6*cff4)
              cff=(1.0_r8+Bio(i,k,iLphy)*cff6)*cff3*cff5*Bio(i,k,iPzoo)
# endif
              Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff)
              GraPL2ZP=cff*Bio(i,k,iLphy)
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraPL2ZP
#endif
#ifdef IRON_LIMIT
!     Grazing Fe terms
              cffFe=GraPL2ZP*Bio(i,k,iFeLp)/MAX(MinVal,Bio(i,k,iLphy))
              Bio(i,k,iFeLp)=Bio(i,k,iFeLp)-cffFe
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+cffFe*FeRR(ng)
#endif
!
!  Predactory Zooplankton grazing on Small Zooplankton, GraZS2ZP.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(ZS2ZPstar(ng)-Bio(i,k,iSzoo)))
              cff5=EXP(-PusaiZS(ng)*Bio(i,k,iLzoo))
              GraZS2ZP=fac7*cff3*cff5*MAX(0.0_r8,cff4)*Bio(i,k,iPzoo)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)-GraZS2ZP
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraZS2ZP
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KZS2ZP(ng)+Bio(i,k,iSzoo)*Bio(i,k,iSzoo))
              cff5=EXP(-PusaiZS(ng)*Bio(i,k,iLzoo))
              cff=fac7*cff3*cff4*cff5*Bio(i,k,iPzoo)*Bio(i,k,iSzoo)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(ZS2ZPstar(ng)-Bio(i,k,iSzoo)))
              cff5=EXP(-PusaiZS(ng)*Bio(i,k,iLzoo))
              cff6=1.0_r8/(fac7*cff4)
              cff=(1.0_r8+Bio(i,k,iSzoo)*cff6)*cff3*cff5*Bio(i,k,iPzoo)
# endif
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)/(1.0_r8+cff)
              GraZS2ZP=cff*Bio(i,k,iSzoo)
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraZS2ZP
#endif
!
!  Predactory Zooplankton grazing on Large Zooplankton, GraZL2ZP.
!
#if defined IVLEV_EXPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(ZL2ZPstar(ng)-Bio(i,k,iLzoo)))
              GraZL2ZP=fac8*cff3*MAX(0.0_r8,cff4)*Bio(i,k,iPzoo)
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)-GraZL2ZP
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraZL2ZP
#else
# ifdef HOLLING_GRAZING
              cff4=1.0_r8/(KZL2ZP(ng)+Bio(i,k,iLzoo)*Bio(i,k,iLzoo))
              cff=fac8*cff3*cff4*Bio(i,k,iPzoo)*Bio(i,k,iLzoo)
# elif defined IVLEV_IMPLICIT
              cff4=1.0_r8-EXP(LamP(ng)*(ZL2ZPstar(ng)-Bio(i,k,iLzoo)))
              cff5=1.0_r8/(fac8*cff4)
              cff=(1.0_r8+Bio(i,k,iLzoo)*cff5)*cff3*Bio(i,k,iPzoo)
# endif
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)/(1.0_r8+cff)
              GraZL2ZP=cff*Bio(i,k,iLzoo)
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)+GraZL2ZP
#endif
!
!  Zooplankton egestion to Particulate Organic Nitrogen (PON) and
!  Particulate Organic Silica (opal).
!
              EgeZS=(1.0_r8-AlphaZS(ng))*(GraPS2ZS+GraPL2ZS)
              EgeZL=(1.0_r8-AlphaZL(ng))*(GraPS2ZL+GraPL2ZL+GraZS2ZL)
              EgeZP=(1.0_r8-AlphaZP(ng))*(GraPL2ZP+GraZS2ZP+GraZL2ZP)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)-EgeZS
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)-EgeZL
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)-EgeZP
              Bio(i,k,iPON_)=Bio(i,k,iPON_)+EgeZS+EgeZL+EgeZP
              Bio(i,k,iopal)=Bio(i,k,iopal)+                            &
     &                       (GraPL2ZS+GraPL2ZL+GraPL2ZP)*LPSiN(i,k)
!
!  Zooplankton excretion to NH4.
!
              ExcZS=(AlphaZS(ng)-BetaZS(ng))*(GraPS2ZS+GraPL2ZS)
              ExcZL=(AlphaZL(ng)-BetaZL(ng))*                           &
     &              (GraPS2ZL+GraPL2ZL+GraZS2ZL)
              ExcZP=(AlphaZP(ng)-BetaZP(ng))*                           &
     &              (GraPL2ZP+GraZS2ZP+GraZL2ZP)
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)-ExcZS
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)-ExcZL
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)-ExcZP
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+ExcZS+ExcZL+ExcZP
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton mortality to particulate organic nitrogen.
!-----------------------------------------------------------------------
!
          fac1=dtdays*MorZS0(ng)
          fac2=dtdays*MorZL0(ng)
          fac3=dtdays*MorZP0(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*Bio(i,k,iSzoo)*EXP(KMorZS(ng)*Bio(i,k,itemp))
              cff2=fac2*Bio(i,k,iLzoo)*EXP(KMorZL(ng)*Bio(i,k,itemp))
              cff3=fac3*Bio(i,k,iPzoo)*EXP(KMorZP(ng)*Bio(i,k,itemp))
              Bio(i,k,iSzoo)=Bio(i,k,iSzoo)/(1.0_r8+cff1)
              Bio(i,k,iLzoo)=Bio(i,k,iLzoo)/(1.0_r8+cff2)
              Bio(i,k,iPzoo)=Bio(i,k,iPzoo)/(1.0_r8+cff3)
              Bio(i,k,iPON_)=Bio(i,k,iPON_)+                            &
     &                       Bio(i,k,iSzoo)*cff1+                       &
     &                       Bio(i,k,iLzoo)*cff2+                       &
     &                       Bio(i,k,iPzoo)*cff3
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Nutrient decomposition.
!-----------------------------------------------------------------------
!
          fac1=dtdays*Nit0(ng)
          fac2=dtdays*VP2N0(ng)
          fac3=dtdays*VP2D0(ng)
          fac4=dtdays*VD2N0(ng)
          fac5=dtdays*VO2S0(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  Nitrification: NH4 to NO3.
!
              cff1=fac1*EXP(KNit(ng)*Bio(i,k,itemp))
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff1)
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iNH4_)*cff1
!
!  Decomposition: PON to NH4.
!
              cff2=fac2*EXP(KP2N(ng)*Bio(i,k,itemp))
              Bio(i,k,iPON_)=Bio(i,k,iPON_)/(1.0_r8+cff2)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+                            &
     &                       Bio(i,k,iPON_)*cff2
!
!  Decomposition: PON to DON.
!
              cff3=fac3*EXP(KP2D(ng)*Bio(i,k,itemp))
              Bio(i,k,iPON_)=Bio(i,k,iPON_)/(1.0_r8+cff3)
              Bio(i,k,iDON_)=Bio(i,k,iDON_)+                            &
     &                       Bio(i,k,iPON_)*cff3
!
!  Decomposition: DON to NH4.
!
              cff4=fac4*EXP(KD2N(ng)*Bio(i,k,itemp))
              Bio(i,k,iDON_)=Bio(i,k,iDON_)/(1.0_r8+cff4)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+                            &
     &                       Bio(i,k,iDON_)*cff4
!
!  Decomposition: Opal to SiOH4.
!
              cff5=fac5*EXP(KO2S(ng)*Bio(i,k,itemp))
              Bio(i,k,iopal)=Bio(i,k,iopal)/(1.0_r8+cff5)
              Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+                            &
     &                       Bio(i,k,iopal)*cff5
            END DO
          END DO
#ifdef NEMURO_SED1
!ENC:  using sediment pools

          DO i=Istr,Iend
!
!  Decomposition: Sediment Opal to SiOH4.
!
            cff6=fac5*EXP(KO2S(ng)*Bio(i,1,itemp))
            OPALsed(i,j)=OPALsed(i,j)/(1.0_r8+cff6)
            Bio(i,1,iSiOH)=Bio(i,1,iSiOH)+                              &
      &                       OPALsed(i,j)*cff6*Hz_inv(i,1)
!
!  Decomposition: Sediment PON to NH4.
!
            cff7=fac2*EXP(KP2N(ng)*Bio(i,1,itemp))
            PONsed(i,j)=PONsed(i,j)/(1.0_r8+cff7)
            Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+                              &
      &                       PONsed(i,j)*cff7*Hz_inv(i,1)
!
!  ENC+CAS:  Our attempt at denitrification
!  This is the NO3 that is lost to N2 due to denitrification
!  in the sediments.  We assume 12% of PONsed remineralization is
!  from denitrification (blame Katja, 2009)
!
            cff8= 0.12_r8*PONsed(i,j)*cff7*5.3_r8
            DENITsed(i,j)=MIN( Bio(i,1,iNO3_)*Hz(i,j,1),cff8)*          &
     &                     (1.0_r8/dtdays)
            Bio(i,1,iNO3_)=MAX(Bio(i,1,iNO3_)-cff8*Hz_inv(i,1),0.0_r8)

          END DO
#endif
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: PON and Opal.
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
                qc(i,k)=Bio(i,k,ibio)
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

#ifdef BIO_SEDIMENT
!
!  Particulate fluxes reaching the seafloor are remineralized and returned
!  to the dissolved nutrient pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nutrient conservation. This may require a time delay
!  remineralization in the future.
!
!  HGA: The original Nemuro model has a restoring upwelling rate (UPW).
!       The code below is an interpretation in terms of the semi-
!       Lagrangian algorithm. What is the correct nutrient path from
!       the benthos to the water column?  NH4 to NO3?
!
!
            IF (ibio.eq.iPON_) THEN
              DO i=Istr,Iend
# ifdef NEMURO_SED1
                cff1=FC(i,0)*6.625_r8*(1.0_r8/dtdays)
                frac_buried(i) = 0.013_r8 + 0.53_r8*cff1**2/            &
     &                     (7.0_r8+cff1)**2
                cff1=FC(i,0)
                PON_burial(i,j)=frac_buried(i)*FC(i,0)*                 &
     &                           (1.0_r8/dtdays)
                PONsed(i,j)=PONsed(i,j)+(1.0_r8-frac_buried(i))*cff1
# else
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iNO3_)=Bio(i,1,iNO3_)+cff1
!ENC: put the sinking PON into the ammonia pool. Do not need this if using
!     the sediment pools
!               Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
!
# endif
              END DO
            ELSE IF (ibio.eq.iopal) THEN
              DO i=Istr,Iend
# ifdef NEMURO_SED1
                cff1=FC(i,0)
                OPAL_burial(i,j)=frac_buried(i)*FC(i,0)*                &
     &                               (1.0_r8/dtdays)
                OPALsed(i,j)=OPALsed(i,j)+(1.0_r8-frac_buried(i))*cff1
# else
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iSiOH)=Bio(i,1,iSiOH)+cff1
# endif
              END DO
            END IF
#endif

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
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO

#ifdef PRIMARY_PROD
        DO i=Istr,Iend
          Bio_NPP(i,j) = 0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio_NPP(i,j) = Bio_NPP(i,j) + Hz(i,j,k)*NPP_slice(i,k)
          END DO
        END DO
#endif

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
