      SUBROUTINE ad_biology (ng,tile)
!
!svn $Id: ad_npzd_iron.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Nutrient-Phytoplankton-Zooplankton-Detritus Model.                  !
!  including Iron Limitation on Phytoplankton Growth.                  !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fiechter, J., A.M. Moore, C.A. Edwards, K.W. Bruland,             !
!      E. Di Lorenzo, C.V.W. Lewis, T.M. Powell, E. Curchitser         !
!      and K. Hedstrom, 2009: Modeling iron limitation of primary      !
!      production in the coastal Gulf of Alaska, Deep Sea Res. II      !
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
      IF (Lbiofile(iADM)) THEN
#else
      IF (Lbiofile(iADM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iADM)=.FALSE.
        BIONAME(iADM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iADM, 15)
#endif
      CALL ad_biology_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#if defined IRON_LIMIT && defined IRON_RELAX
     &                      GRID(ng) % h,                               &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % ad_z_r,                          &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % ad_z_w,                          &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % ad_srflx,                      &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t)

#ifdef PROFILE
      CALL wclock_off (ng, iADM, 15)
#endif
      RETURN
      END SUBROUTINE ad_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE ad_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
#if defined IRON_LIMIT && defined IRON_RELAX
     &                            h,                                    &
#endif
     &                            Hz, ad_Hz, z_r, ad_z_r, z_w, ad_z_w,  &
     &                            srflx, ad_srflx,                      &
     &                            t, ad_t)
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
#if defined IRON_LIMIT && defined IRON_RELAX
      real(r8), intent(in) :: h(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_srflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
#if defined IRON_LIMIT && defined IRON_RELAX
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: ad_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: ad_srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 2

      integer :: Iter, i, ibio, isink, itime, itrc, iTrcMax, j, k, ks
      integer :: Iteradj, kk

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: MinVal = 1.0e-6_r8

      real(r8) :: Att, ExpAtt, Itop, PAR, PAR1
      real(r8) :: ad_Att, ad_ExpAtt, ad_Itop, ad_PAR
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, dtdays
      real(r8) :: ad_cff, ad_cff1, ad_cff4, ad_cff5, ad_cff6
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: ad_cffL, ad_cffR, ad_cu, ad_dltL, ad_dltR
      real(r8) :: fac, fac1, fac2
      real(r8) :: ad_fac, ad_fac1, ad_fac2
      real(r8) :: adfac, adfac1, adfac2, adfac3
#ifdef IRON_LIMIT
      real(r8) :: Nlimit, FNlim
      real(r8) :: ad_Nlimit, ad_FNlim
      real(r8) :: FNratio, FCratio, FCratioE, Flimit
      real(r8) :: ad_FNratio, ad_FCratio, ad_FCratioE, ad_Flimit
      real(r8) :: FeC2FeN, FeN2FeC
# ifdef IRON_RELAX
      real(r8) :: FeNudgCoef
# endif
#endif
      real(r8), dimension(Nsink) :: Wbio
      real(r8), dimension(Nsink) :: ad_Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: ad_PARsur

      real(r8), dimension(NT(ng),2) :: BioTrc
      real(r8), dimension(NT(ng),2) :: BioTrc1
      real(r8), dimension(NT(ng),2) :: ad_BioTrc
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio1
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio2
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Light
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Light
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_qc

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

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Set nudging coefficient for dissolved iron over the shelf.
!
      FeNudgCoef=dt(ng)/(FeNudgTime(ng)*86400.0_r8)
#endif
#ifdef IRON_LIMIT
!
!  Set Fe:N and Fe:C conversion ratio and its inverse.
!
          FeN2FeC=(16.0_r8/106.0_r8)*1.0E3_r8
          FeC2FeN=(106.0_r8/16.0_r8)*1.0E-3_r8
#endif
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt                 ! Phytoplankton
      idsink(2)=iSdet                 ! Small detritus
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! Phytoplankton
      Wbio(2)=wDet(ng)                ! Small detritus
!
      ad_Wbio(1)=0.0_r8
      ad_Wbio(2)=0.0_r8
!
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
        ad_PAR=0.0_r8
        ad_Att=0.0_r8
        ad_ExpAtt=0.0_r8
        ad_Itop=0.0_r8
        ad_dltL=0.0_r8
        ad_dltR=0.0_r8
        ad_cu=0.0_r8
        ad_cff=0.0_r8
        ad_cff1=0.0_r8
        ad_cff4=0.0_r8
        ad_cff5=0.0_r8
        ad_cff6=0.0_r8
        ad_fac=0.0_r8
        ad_fac1=0.0_r8
        ad_fac2=0.0_r8
        ad_cffL=0.0_r8
        ad_cffR=0.0_r8
        adfac=0.0_r8
        adfac1=0.0_r8
        adfac2=0.0_r8
        adfac3=0.0_r8
#ifdef IRON_LIMIT
        ad_FNratio=0.0_r8
        ad_FCratio=0.0_r8
        ad_FCratioE=0.0_r8
        ad_Flimit=0.0_r8
        ad_FNlim=0.0_r8
        ad_Nlimit=0.0_r8
#endif
!
        DO k=1,N(ng)
          DO i=IminS,ImaxS
            ad_Hz_inv(i,k)=0.0_r8
            ad_Hz_inv2(i,k)=0.0_r8
            ad_Hz_inv3(i,k)=0.0_r8
            ad_WL(i,k)=0.0_r8
            ad_WR(i,k)=0.0_r8
            ad_bL(i,k)=0.0_r8
            ad_bR(i,k)=0.0_r8
            ad_qc(i,k)=0.0_r8
            ad_Light(i,k)=0.0_r8
          END DO
        END DO
        DO itrc=1,NBT
          ibio=idbio(itrc)
          ad_BioTrc(ibio,1)=0.0_r8
          ad_BioTrc(ibio,2)=0.0_r8
        END DO
        DO i=IminS,ImaxS
          ad_PARsur(i)=0.0_r8
        END DO
        DO k=0,N(ng)
          DO i=IminS,ImaxS
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
!
!  Clear ad_Bio and Bio arrays.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,ibio)=0.0_r8
              Bio1(i,k,ibio)=0.0_r8
              Bio2(i,k,ibio)=0.0_r8
              ad_Bio(i,k,ibio)=0.0_r8
              ad_Bio_old(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
!
!  Compute inverse thickness to avoid repeated divisions.
!
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
!  Compute the required basic state arrays.
!
!  Restrict biological tracer to be positive definite. If a negative
!  concentration is detected, nitrogen is drawn from the most abundant
!  pool to supplement the negative pools to a lower limit of MinVal
!  which is set to 1E-6 above.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
!  NOTE: In the following code, t(:,:,:,nnew,:) should be in units of
!        tracer times depth. However the basic state (nstp and nnew
!        indices) that is read from the forward file is in units of
!        tracer. Since BioTrc(ibio,nnew) is in tracer units, we simply
!        use t instead of t*Hz_inv.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
!>            BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>
              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>            BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
!>
              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)
            END DO
!
!  Impose positive definite concentrations.
!
            cff2=0.0_r8
            DO itime=1,2
              cff1=0.0_r8
              iTrcMax=idbio(1)
#ifdef IRON_LIMIT
              DO itrc=1,NBT-2
#else
              DO itrc=1,NBT
#endif
                ibio=idbio(itrc)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                  iTrcMax=ibio
                END IF
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
              IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
              END IF
#ifdef IRON_LIMIT
              DO itrc=NBT-1,NBT
                ibio=idbio(itrc)
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
#endif
            END DO
!
!  Load biological tracers into local arrays.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
              Bio_old(i,k,ibio)=BioTrc(ibio,nstp)
              Bio(i,k,ibio)=BioTrc(ibio,nstp)
            END DO

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
            IF (h(i,j).le.FeHmin(ng)) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
            END IF
#endif
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
#ifdef CONST_PAR
!
!  Specify constant surface irradiance a la Powell and Spitz.
!
          PARsur(i)=158.075_r8
#else
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#endif
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
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                Att=(AttSW(ng)+AttPhy(ng)*Bio(i,k,iPhyt))*              &
     &              (z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                Light(i,k)=0.0_r8
              END DO
            END IF
          END DO
!
!  Phytoplankton photosynthetic growth and nitrate uptake (Vm_NO3 rate).
!  The Michaelis-Menten curve is used to describe the change in uptake
!  rate as a function of nitrate concentration. Here, PhyIS is the
!  initial slope of the P-I curve and K_NO3 is the half saturation of
!  phytoplankton nitrate uptake.
#ifdef IRON_LIMIT
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
#endif
!
          cff1=dtdays*Vm_NO3(ng)*PhyIS(ng)
          cff2=Vm_NO3(ng)*Vm_NO3(ng)
          cff3=PhyIS(ng)*PhyIS(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef IRON_LIMIT
              FNratio=Bio(i,k,iFphy)/MAX(MinVal,Bio(i,k,iPhyt))
              FCratio=FNratio*FeN2FeC
              FCratioE=B_Fe(ng)*Bio(i,k,iFdis)**A_Fe(ng)
              Flimit=FCratio*FCratio/                                   &
     &               (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))

              Nlimit=1.0_r8/(K_NO3(ng)+Bio(i,k,iNO3_))
              FNlim=MIN(1.0_r8,Flimit/(Bio(i,k,iNO3_)*Nlimit))
#endif
              cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
              cff=Bio(i,k,iPhyt)*                                       &
#ifdef IRON_LIMIT
     &            cff1*cff4*Light(i,k)*FNlim*Nlimit
#else

     &            cff1*cff4*Light(i,k)/                                 &
     &            (K_NO3(ng)+Bio(i,k,iNO3_))
#endif

              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                            &
     &                       Bio(i,k,iNO3_)*cff
#ifdef IRON_LIMIT
!
!  Iron uptake proportional to growth.
!
              fac=cff*Bio(i,k,iNO3_)*FNratio/                           &
     &            MAX(MinVal,Bio(i,k,iFdis))
              Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+fac)
              Bio(i,k,iFphy)=Bio(i,k,iFphy)+                            &
     &                       Bio(i,k,iFdis)*fac
!
!  Iron uptake to reach appropriate Fe:C ratio.
!
              cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
              cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
              IF (cff6.ge.0.0_r8) THEN
                cff=cff6/MAX(MinVal,Bio(i,k,iFdis))
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*cff
              ELSE
                cff=-cff6/MAX(MinVal,Bio(i,k,iFphy))
                Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*cff
              END IF
#endif
            END DO
          END DO
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
#ifdef IRON_LIMIT
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
#endif
!
          cff1=dtdays*ZooGR(ng)
          cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,iZoop)*                                       &
     &            cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio(i,k,iPhyt)))/         &
     &            Bio(i,k,iPhyt)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/(1.0_r8+cff)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+                            &
     &                       Bio(i,k,iPhyt)*cff2*cff
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iPhyt)*ZooEEN(ng)*cff
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iPhyt)*ZooEED(ng)*cff
#ifdef IRON_LIMIT
              Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       Bio(i,k,iFphy)*cff*FeRR(ng)
#endif
            END DO
          END DO
!
!  Phytoplankton mortality to nutrients (PhyMRNro rate), detritus
!  (PhyMRD rate), and if applicable dissolved iron (FeRR rate).
!
          cff3=dtdays*PhyMRD(ng)
          cff2=dtdays*PhyMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iPhyt)*cff2
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iPhyt)*cff3
#ifdef IRON_LIMIT
              Bio(i,k,iFphy)=Bio(i,k,iFphy)*cff1
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
#endif
            END DO
          END DO
!
!  Zooplankton mortality to nutrients (ZooMRN rate) and Detritus
!  (ZooMRD rate).
!
          cff3=dtdays*ZooMRD(ng)
          cff2=dtdays*ZooMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iZoop)*cff2
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iZoop)*cff3
            END DO
          END DO
!
!  Detritus breakdown to nutrients: remineralization (DetRR rate).
!
          cff2=dtdays*DetRR(ng)
          cff1=1.0_r8/(1.0_r8+cff2)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iSDet)*cff2
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Phytoplankton and Detritus
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
              bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
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
!>            tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
!>   &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)
!>
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+cff*ad_t(i,j,k,nnew,ibio)
              ad_cff=ad_cff+Hz(i,j,k)*ad_t(i,j,k,nnew,ibio)
!>            tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
!>
              ad_Bio_old(i,k,ibio)=ad_Bio_old(i,k,ibio)-ad_cff
              ad_Bio(i,k,ibio)=ad_Bio(i,k,ibio)+ad_cff
              ad_cff=0.0_r8
            END DO
          END DO
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
        ITER_LOOP1: DO Iter=BioIter(ng),1,-1
!
!-----------------------------------------------------------------------
!  Adjoint vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
!  Compute appropriate basic state arrays III.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
!  NOTE: In the following code, t(:,:,:,nnew,:) should be in units of
!        tracer times depth. However the basic state (nstp and nnew
!        indices) that is read from the forward file is in units of
!        tracer. Since BioTrc(ibio,nnew) is in tracer units, we simply
!        use t instead of t*Hz_inv.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
!>              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>
                BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
!>
                BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)
              END DO
!
!  Impose positive definite concentrations.
!
              cff2=0.0_r8
              DO itime=1,2
                cff1=0.0_r8
                iTrcMax=idbio(1)
#ifdef IRON_LIMIT
                DO itrc=1,NBT-2
#else
                DO itrc=1,NBT
#endif
                  ibio=idbio(itrc)
                  cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                  IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                    iTrcMax=ibio
                  END IF
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
                IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                  BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
                END IF
#ifdef IRON_LIMIT
                DO itrc=NBT-1,NBT
                  ibio=idbio(itrc)
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
#endif
              END DO
!
!  Load biological tracers into local arrays.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
                Bio_old(i,k,ibio)=BioTrc(ibio,nstp)
                Bio(i,k,ibio)=BioTrc(ibio,nstp)
              END DO

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
              IF (h(i,j).le.FeHmin(ng)) THEN
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
              END IF
#endif
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
#ifdef CONST_PAR
!
!  Specify constant surface irradiance a la Powell and Spitz.
!
            PARsur(i)=158.075_r8
#else
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#endif
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Compute light attenuation as function of depth.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              IF (PARsur(i).gt.0.0_r8) THEN            ! day time
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                  Att=(AttSW(ng)+AttPhy(ng)*Bio(i,k,iPhyt))*            &
     &                (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att  ! average at cell center
                  Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
              ELSE                                     ! night time
                DO k=1,N(ng)
                  Light(i,k)=0.0_r8
                END DO
              END IF
            END DO
!
!  Phytoplankton photosynthetic growth and nitrate uptake (Vm_NO3 rate).
!  The Michaelis-Menten curve is used to describe the change in uptake
!  rate as a function of nitrate concentration. Here, PhyIS is the
!  initial slope of the P-I curve and K_NO3 is the half saturation of
!  phytoplankton nitrate uptake.
#ifdef IRON_LIMIT
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
#endif
!
            cff1=dtdays*Vm_NO3(ng)*PhyIS(ng)
            cff2=Vm_NO3(ng)*Vm_NO3(ng)
            cff3=PhyIS(ng)*PhyIS(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
#ifdef IRON_LIMIT
!
!  Calculate growth reduction factor due to iron limitation.
!
                FNratio=Bio(i,k,iFphy)/MAX(MinVal,Bio(i,k,iPhyt))
                FCratio=FNratio*FeN2FeC
                FCratioE=B_Fe(ng)*Bio(i,k,iFdis)**A_Fe(ng)
                Flimit=FCratio*FCratio/                                 &
     &                 (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))

                Nlimit=1.0_r8/(K_NO3(ng)+Bio(i,k,iNO3_))
                FNlim=MIN(1.0_r8,Flimit/(Bio(i,k,iNO3_)*Nlimit))
#endif
                cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
                cff=Bio(i,k,iPhyt)*                                     &
#ifdef IRON_LIMIT
     &              cff1*cff4*Light(i,k)*FNlim*Nlimit
#else
     &              cff1*cff4*Light(i,k)/                               &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
#endif
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff
#ifdef IRON_LIMIT
!
!  Iron uptake proportional to growth.
!
                fac=cff*Bio(i,k,iNO3_)*FNratio/                         &
     &              MAX(MinVal,Bio(i,k,iFdis))
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+fac)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*fac


! Iron Uptake proportional to growth
                cff=cff*Bio(i,k,iNO3_)*                                 &
     &              FNratio/max(MinVal,Bio(i,k,iFdis))
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*cff
!
!  Iron uptake to reach appropriate Fe:C ratio.
!
                cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
                cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
                IF (cff6.ge.0.0_r8) THEN
                  cff=cff6/MAX(MinVal,Bio(i,k,iFdis))
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)+                        &
     &                           Bio(i,k,iFdis)*cff
                ELSE
                  cff=-cff6/MAX(MinVal,Bio(i,k,iFphy))
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)+                        &
     &                           Bio(i,k,iFphy)*cff
                END IF
#endif
              END DO
            END DO
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
#ifdef IRON_LIMIT
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
#endif
!
            cff1=dtdays*ZooGR(ng)
            cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iZoop)*                                     &
     &              cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio(i,k,iPhyt)))/       &
     &              Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/(1.0_r8+cff)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         Bio(i,k,iPhyt)*cff2*cff
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*ZooEEN(ng)*cff
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*ZooEED(ng)*cff
#ifdef IRON_LIMIT
                Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*cff*FeRR(ng)
#endif
              END DO
            END DO
!
!  Phytoplankton mortality to nutrients (PhyMRNro rate), detritus
!  (PhyMRD rate), and if applicable dissolved iron (FeRR rate).
!
            cff3=dtdays*PhyMRD(ng)
            cff2=dtdays*PhyMRN(ng)
            cff1=1.0_r8/(1.0_r8+cff2+cff3)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff1
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*cff2
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*cff3
#ifdef IRON_LIMIT
                Bio(i,k,iFphy)=Bio(i,k,iFphy)*cff1
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
#endif
              END DO
            END DO
!
!  Zooplankton mortality to nutrients (ZooMRN rate) and Detritus
!  (ZooMRD rate).
!
            cff3=dtdays*ZooMRD(ng)
            cff2=dtdays*ZooMRN(ng)
            cff1=1.0_r8/(1.0_r8+cff2+cff3)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iZoop)*cff2
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iZoop)*cff3
              END DO
            END DO
!
!  Detritus breakdown to nutrients: remineralization (DetRR rate).
!
            cff2=dtdays*DetRR(ng)
            cff1=1.0_r8/(1.0_r8+cff2)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff1
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iSDet)*cff2
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Phytoplankton and Detritus
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
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
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
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
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays III.
!
          SINK_LOOP1: DO isink=1,Nsink
            ibio=idsink(isink)

!
!  Compute required flux arrays.
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
              bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
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
!>              tl_Bio(i,k,ibio)=tl_qc(i,k)+                            &
!>   &                           (tl_FC(i,k)-tl_FC(i,k-1))*Hz_inv(i,k)+ &
!>   &                           (FC(i,k)-FC(i,k-1))*tl_Hz_inv(i,k)
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_Bio(i,k,ibio)
                ad_FC(i,k)=ad_FC(i,k)+Hz_inv(i,k)*ad_Bio(i,k,ibio)
                ad_FC(i,k-1)=ad_FC(i,k-1)-Hz_inv(i,k)*ad_Bio(i,k,ibio)
                ad_Hz_inv(i,k)=ad_Hz_inv(i,k)+                          &
     &                         (FC(i,k)-FC(i,k-1))*ad_Bio(i,k,ibio)
                ad_Bio(i,k,ibio)=0.0_r8
              END DO
            END DO
!
!  Adjoint of final computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
!>              tl_FC(i,k-1)=tl_FC(i,k-1)+                              &
!>   &                       (tl_Hz(i,j,ks)*cu+Hz(i,j,ks)*tl_cu)*       &
!>   &                       (bL(i,ks)+                                 &
!>   &                        cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-           &
!>   &                            (1.5_r8-cu)*                          &
!>   &                            (bR(i,ks)+bL(i,ks)-                   &
!>   &                             2.0_r8*qc(i,ks))))+                  &
!>   &                       Hz(i,j,ks)*cu*                             &
!>   &                       (tl_bL(i,ks)+                              &
!>   &                        tl_cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-        &
!>   &                               (1.5_r8-cu)*                       &
!>   &                               (bR(i,ks)+bL(i,ks)-                &
!>   &                               2.0_r8*qc(i,ks)))+                 &
!>   &                        cu*(0.5_r8*(tl_bR(i,ks)-tl_bL(i,ks))+     &
!>   &                            tl_cu*                                &
!>   &                            (bR(i,ks)+bL(i,ks)-2.0_r8*qc(i,ks))-  &
!>   &                            (1.5_r8-cu)*                          &
!>   &                            (tl_bR(i,ks)+tl_bL(i,ks)-             &
!>   &                             2.0_r8*tl_qc(i,ks))))
!>
                adfac=(bL(i,ks)+                                        &
     &                 cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-                  &
     &                     (1.5_r8-cu)*                                 &
     &                     (bR(i,ks)+bL(i,ks)-                          &
     &                      2.0_r8*qc(i,ks))))*ad_FC(i,k-1)
                adfac1=Hz(i,j,ks)*cu*ad_FC(i,k-1)
                adfac2=adfac1*cu
                adfac3=adfac2*(1.5_r8-cu)
                ad_Hz(i,j,ks)=ad_Hz(i,j,ks)+cu*adfac
                ad_cu=ad_cu+Hz(i,j,ks)*adfac
                ad_bL(i,ks)=ad_bL(i,ks)+adfac1
                ad_cu=ad_cu+                                            &
     &                adfac1*(0.5_r8*(bR(i,ks)-bL(i,ks))-               &
     &                        (1.5_r8-cu)*                              &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks)))+                      &
     &                adfac2*(bR(i,ks)+bL(i,ks)-2.0_r8*qc(i,ks))
                ad_bR(i,ks)=ad_bR(i,ks)+0.5_r8*adfac2-adfac3
                ad_bL(i,ks)=ad_bL(i,ks)-0.5_r8*adfac2-adfac3
                ad_qc(i,ks)=ad_qc(i,ks)+2.0_r8*adfac3
!>              tl_cu=(0.5_r8+SIGN(0.5_r8,                              &
!>   &                             (1.0_r8-(WL(i,k)-z_w(i,j,ks-1))*     &
!>   &                             Hz_inv(i,ks))))*                     &
!>   &                ((tl_WL(i,k)-tl_z_w(i,j,ks-1))*Hz_inv(i,ks)+      &
!>   &                 (WL(i,k)-z_w(i,j,ks-1))*tl_Hz_inv(i,ks))
!>
                adfac=(0.5_r8+SIGN(0.5_r8,                              &
     &                             (1.0_r8-(WL(i,k)-z_w(i,j,ks-1))*     &
     &                             Hz_inv(i,ks))))*ad_cu
                adfac1=adfac*Hz_inv(i,ks)
                ad_WL(i,k)=ad_WL(i,k)+adfac1
                ad_z_w(i,j,ks-1)=ad_z_w(i,j,ks-1)-adfac1
                ad_Hz_inv(i,ks)=ad_Hz_inv(i,ks)+                        &
     &                          (WL(i,k)-z_w(i,j,ks-1))*adfac
                ad_cu=0.0_r8
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
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
!>                  tl_FC(i,k-1)=tl_FC(i,k-1)+tl_WR(i,ks)
!>
                    ad_WR(i,ks)=ad_WR(i,ks)+ad_FC(i,k-1)
                  END IF
                END DO
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_WR(i,k)=tl_Hz(i,j,k)*qc(i,k)+Hz(i,j,k)*tl_qc(i,k)
!>
                ad_Hz(i,j,k)=ad_Hz(i,j,k)+qc(i,k)*ad_WR(i,k)
                ad_qc(i,k)=ad_qc(i,k)+Hz(i,j,k)*ad_WR(i,k)
                ad_WR(i,k)=0.0_r8
!>              tl_WL(i,k)=tl_z_w(i,j,k-1)+tl_cff
!>
                ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)+ad_WL(i,k)
                ad_cff=ad_cff+ad_WL(i,k)
                ad_WL(i,k)=0.0_r8
!>              tl_FC(i,k-1)=0.0_r8
!>
                ad_FC(i,k-1)=0.0_r8
              END DO
            END DO
!>          tl_cff=dtdays*SIGN(1.0_r8,Wbio(isink))*tl_Wbio(isink)
!>
            ad_Wbio(isink)=ad_Wbio(isink)+                              &
     &                     dtdays*SIGN(1.0_r8,Wbio(isink))*ad_cff
            ad_cff=0.0_r8
!
!  Compute appropriate values of bR and bL.
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
              bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
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
!>              tl_bL(i,k)=tl_qc(i,k)-tl_dltL
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_bL(i,k)
                ad_dltL=ad_dltL-ad_bL(i,k)
                ad_bL(i,k)=0.0_r8
!>              tl_bR(i,k)=tl_qc(i,k)+tl_dltR
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_bR(i,k)
                ad_dltR=ad_dltR+ad_bR(i,k)
                ad_bR(i,k)=0.0_r8
                IF ((dltR*dltL).lt.0.0_r8) THEN
!>                tl_dltR=0.0_r8
!>
                  ad_dltR=0.0_r8
!>                tl_dltL=0.0_r8
!>
                  ad_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
!>                tl_dltR=tl_cffL
!>
                  ad_cffL=ad_cffL+ad_dltR
                  ad_dltR=0.0_r8
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
!>                tl_dltL=tl_cffR
!>
                  ad_cffR=ad_cffR+ad_dltL
                  ad_dltL=0.0_r8
                END IF
!>              tl_cffL=2.0_r8*tl_dltL
!>
                ad_dltL=ad_dltL+2.0_r8*ad_cffL
                ad_cffL=0.0_r8
!>              tl_cffR=2.0_r8*tl_dltR
!>
                ad_dltR=ad_dltR+2.0_r8*ad_cffR
                ad_cffR=0.0_r8
!>              tl_dltL=tl_qc(i,k)-tl_bL(i,k)
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_dltL
                ad_bL(i,k)=ad_bL(i,k)-ad_dltL
                ad_dltL=0.0_r8
!>              tl_dltR=tl_bR(i,k)-tl_qc(i,k)
!>
                ad_bR(i,k)=ad_bR(i,k)+ad_dltR
                ad_qc(i,k)=ad_qc(i,k)-ad_dltR
                ad_dltR=0.0_r8
              END DO
            END DO
            DO i=Istr,Iend
#if defined LINEAR_CONTINUATION
!>            tl_bR(i,1)=tl_bL(i,2)
!>
              ad_bL(i,2)=ad_bL(i,2)+ad_bR(i,1)
              ad_bR(i,1)=0.0_r8
!>            tl_bL(i,1)=2.0_r8*tl_qc(i,1)-tl_bR(i,1)
!>
              ad_qc(i,1)=ad_qc(i,1)+2.0_r8*ad_bL(i,1)
              ad_bR(i,1)=ad_bR(i,1)-ad_bL(i,1)
              ad_bL(i,1)=0.0_r8
#elif defined NEUMANN
!>            tl_bR(i,1)=tl_bL(i,2)
!>
              ad_bL(i,2)=ad_bL(i,2)+ad_bR(i,1)
              ad_bR(i,1)=0.0_r8
!>            tl_bL(i,1)=1.5_r8*tl_qc(i,1)-0.5_r8*tl_bR(i,1)
!>
              ad_qc(i,1)=ad_qc(i,1)+1.5_r8*ad_bL(i,1)
              ad_bR(i,1)=ad_bR(i,1)-0.5_r8*ad_bL(i,1)
              ad_bL(i,1)=0.0_r8
#else
!>            tl_bL(i,2)=tl_qc(i,1)         ! bottom grid boxes are
!>            tl_bR(i,1)=tl_qc(i,1)         ! re-assumed to be
!>            tl_bL(i,1)=tl_qc(i,1)         ! piecewise constant.
!>
              ad_qc(i,1)=ad_qc(i,1)+ad_bL(i,1)+                         &
     &                              ad_bR(i,1)+                         &
     &                              ad_bL(i,2)
              ad_bL(i,1)=0.0_r8
              ad_bR(i,1)=0.0_r8
              ad_bL(i,2)=0.0_r8
#endif
#if defined LINEAR_CONTINUATION
!>            tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
!>
              ad_bR(i,N(ng)-1)=ad_bR(i,N(ng)-1)+ad_bL(i,N(ng))
              ad_bL(i,N(ng))=0.0_r8
!>            tl_bR(i,N(ng))=2.0_r8*tl_qc(i,N(ng))-tl_bL(i,N(ng))
!>
              ad_qc(i,N(ng))=ad_qc(i,N(ng))+2.0_r8*ad_bR(i,N(ng))
              ad_bL(i,N(ng))=ad_bL(i,N(ng))-ad_bR(i,N(ng))
              ad_bR(i,N(ng))=0.0_r8
#elif defined NEUMANN
!>            tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
!>
              ad_bR(i,N(ng)-1)=ad_bR(i,N(ng)-1)+ad_bL(i,N(ng))
              ad_bL(i,N(ng))=0.0_r8
!>            tl_bR(i,N(ng))=1.5_r8*tl_qc(i,N(ng))-0.5_r8*tl_bL(i,N(ng))
!>
              ad_qc(i,N(ng))=ad_qc(i,N(ng))+1.5_r8*ad_bR(i,N(ng))
              ad_bL(i,N(ng))=ad_bL(i,N(ng))-0.5_r8*ad_bR(i,N(ng))
              ad_bR(i,N(ng))=0.0_r8
#else
!>            tl_bR(i,N(ng))=tl_qc(i,N(ng)) ! default strictly monotonic
!>            tl_bL(i,N(ng))=tl_qc(i,N(ng)) ! conditions
!>            tl_bR(i,N(ng)-1)=tl_qc(i,N(ng))
!>
              ad_qc(i,N(ng))=ad_qc(i,N(ng))+ad_bR(i,N(ng)-1)+           &
     &                                      ad_bL(i,N(ng))+             &
     &                                      ad_bR(i,N(ng))
              ad_bR(i,N(ng)-1)=0.0_r8
              ad_bL(i,N(ng))=0.0_r8
              ad_bR(i,N(ng))=0.0_r8
#endif
            END DO
!
!  Compute  WR and WL arrays appropriate for this part of the code.
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
                bR1(i,k)=bR(i,k)
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL1(i,k+1)=bL(i,k+1)
                bL(i,k+1)=bR(i,k)
!>              tl_bL(i,k+1)=tl_bR(i,k)
!>
                ad_bR(i,k)=ad_bR(i,k)+ad_bL(i,k+1)
                ad_bL(i,k+1)=0.0_r8
!>              tl_bR(i,k)=(tl_dltR*bR1(i,k)+dltR*tl_bR(i,k)+           &
!>   &                      tl_dltL*bL1(i,k+1)+dltL*tl_bL(i,k+1))/      &
!>   &                      (dltR+dltL)-                                &
!>   &                      (tl_dltR+tl_dltL)*bR(i,k)/(dltR+dltL)
!>
                adfac=ad_bR(i,k)/(dltR+dltL)
                adfac1=ad_bR(i,k)*bR(i,k)/(dltR+dltL)
                ad_dltR=ad_dltR+adfac*bR1(i,k)
                ad_dltL=ad_dltL+adfac*bL1(i,k+1)
                ad_bL(i,k+1)=ad_bL(i,k+1)+dltL*adfac
                ad_dltR=ad_dltR-adfac1
                ad_dltL=ad_dltL-adfac1
                ad_bR(i,k)=dltR*adfac
!>              tl_dltR=(0.5_r8-SIGN(0.5_r8,cff-WR(i,k+1)))*            &
!>   &                  tl_WR(i,k+1)
!>
                ad_WR(i,k+1)=ad_WR(i,k+1)+                              &
     &                       (0.5_r8-SIGN(0.5_r8,cff-WR(i,k+1)))*       &
     &                       ad_dltR
                ad_dltR=0.0_r8
!>              tl_dltL=(0.5_r8-SIGN(0.5_r8,cff-WL(i,k  )))*            &
!>   &                  tl_WL(i,k  )
!>
                ad_WL(i,k  )=ad_WL(i,k  )+                              &
     &                       (0.5_r8-SIGN(0.5_r8,cff-WL(i,k  )))*       &
     &                       ad_dltL
                ad_dltL=0.0_r8
              END DO
            END DO

            DO k=2,N(ng)-1
              DO i=Istr,Iend
!
!  Compute appropriate dltL and dltr.
!
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
!>              tl_WL(i,k)=2.0_r8*(dltR-2.0_r8*dltL)*                   &
!>   &                            (tl_dltR-2.0_r8*tl_dltL)
!>
                adfac=ad_WL(i,k)*2.0_r8*(dltR-2.0_r8*dltL)
                ad_dltR=ad_dltR+adfac
                ad_dltL=ad_dltL-2.0_r8*adfac
                ad_WL(i,k)=0.0_r8

!>              tl_WR(i,k)=2.0_r8*(2.0_r8*dltR-dltL)*                   &
!>   &                            (2.0_r8*tl_dltR-tl_dltL)
!>
                adfac=ad_WR(i,k)*2.0_r8*(2.0_r8*dltR-dltL)
                ad_dltR=ad_dltR+2.0_r8*adfac
                ad_dltL=ad_dltL-adfac
                ad_WR(i,k)=0.0_r8
!>              tl_bL(i,k)=tl_qc(i,k)-tl_dltL
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_bL(i,k)
                ad_dltL=ad_dltL-ad_bL(i,k)
                ad_bL(i,k)=0.0_r8
!>              tl_bR(i,k)=tl_qc(i,k)+tl_dltR
!>
                ad_qc(i,k)=ad_qc(i,k)+ad_bR(i,k)
                ad_dltR=ad_dltR+ad_bR(i,k)
                ad_bR(i,k)=0.0_r8

!
!  Compute appropriate dltL and dltr.
!
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

                cff=(dltR-dltL)*Hz_inv3(i,k)
!>              tl_dltL=tl_dltL+tl_cff*Hz(i,j,k-1)+cff*tl_Hz(i,j,k-1)
!>
                ad_cff=ad_cff+ad_dltL*Hz(i,j,k-1)
                ad_Hz(i,j,k-1)=ad_Hz(i,j,k-1)+cff*ad_dltL
!>              tl_dltR=tl_dltR-tl_cff*Hz(i,j,k+1)-cff*tl_Hz(i,j,k+1)
!>
                ad_cff=ad_cff-ad_dltR*Hz(i,j,k+1)
                ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)-cff*ad_dltR
!>              tl_cff=(tl_dltR-tl_dltL)*Hz_inv3(i,k)+                  &
!>   &                 (dltR-dltL)*tl_Hz_inv3(i,k)
!>
                adfac=ad_cff*Hz_inv3(i,k)
                ad_dltR=ad_dltR+adfac
                ad_dltL=ad_dltL-adfac
                ad_Hz_inv3(i,k)=ad_Hz_inv3(i,k)+(dltR-dltL)*ad_cff
                ad_cff=0.0_r8
!
!  Compute appropriate dltL and dltr.
!
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)

                IF ((dltR*dltL).le.0.0_r8) THEN
!>                tl_dltR=0.0_r8
!>
                  ad_dltR=0.0_r8
!>                tl_dltL=0.0_r8
!>
                  ad_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
!>                tl_dltR=tl_cffL
!>
                  ad_cffL=ad_cffL+ad_dltR
                  ad_dltR=0.0_r8
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
!>                tl_dltL=tl_cffR
!>
                  ad_cffR=ad_cffR+ad_dltL
                  ad_dltL=0.0_r8
                END IF
!>              tl_cffL=tl_cff*FC(i,k-1)+cff*tl_FC(i,k-1)
!>
                ad_cff=ad_cff+ad_cffL*FC(i,k-1)
                ad_FC(i,k-1)=ad_FC(i,k-1)+cff*ad_cffL
                ad_cffL=0.0_r8
!>              tl_cffR=tl_cff*FC(i,k)+cff*tl_FC(i,k)
!>
                ad_cff=ad_cff+ad_cffR*FC(i,k)
                ad_FC(i,k)=ad_FC(i,k)+cff*ad_cffR
                ad_cffR=0.0_r8
!>              tl_cff=tl_Hz(i,j,k-1)+2.0_r8*tl_Hz(i,j,k)+tl_Hz(i,j,k+1)
!>
                ad_Hz(i,j,k-1)=ad_Hz(i,j,k-1)+ad_cff
                ad_Hz(i,j,k)=ad_Hz(i,j,k)+2.0_r8*ad_cff
                ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)+ad_cff
                ad_cff=0.0_r8
!>              tl_dltL=tl_Hz(i,j,k)*FC(i,k-1)+Hz(i,j,k)*tl_FC(i,k-1)
!>
                ad_Hz(i,j,k)=ad_Hz(i,j,k)+ad_dltL*FC(i,k-1)
                ad_FC(i,k-1)=ad_FC(i,k-1)+ad_dltL*Hz(i,j,k)
                ad_dltL=0.0_r8
!>              tl_dltR=tl_Hz(i,j,k)*FC(i,k)+Hz(i,j,k)*tl_FC(i,k)
!>
                ad_Hz(i,j,k)=ad_Hz(i,j,k)+ad_dltR*FC(i,k)
                ad_FC(i,k)=ad_FC(i,k)+ad_dltR*Hz(i,j,k)
                ad_dltR=0.0_r8
              END DO
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
!>              tl_FC(i,k)=(tl_qc(i,k+1)-tl_qc(i,k))*Hz_inv2(i,k)+      &
!>   &                     (qc(i,k+1)-qc(i,k))*tl_Hz_inv2(i,k)
!>
                adfac=ad_FC(i,k)*Hz_inv2(i,k)
                ad_qc(i,k+1)=ad_qc(i,k+1)+adfac
                ad_qc(i,k)=ad_qc(i,k)-adfac
                ad_Hz_inv2(i,k)=ad_Hz_inv2(i,k)+(qc(i,k+1)-qc(i,k))*    &
     &                          ad_FC(i,k)
                ad_FC(i,k)=0.0_r8
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_qc(i,k)=tl_Bio(i,k,ibio)
!>
                ad_Bio(i,k,ibio)=ad_Bio(i,k,ibio)+ad_qc(i,k)
                ad_qc(i,k)=0.0_r8
              END DO
            END DO

          END DO SINK_LOOP1
!
!  Compute appropriate basic state arrays II.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
!  NOTE: In the following code, t(:,:,:,nnew,:) should be in units of
!        tracer times depth. However the basic state (nstp and nnew
!        indices) that is read from the forward file is in units of
!        tracer. Since BioTrc(ibio,nnew) is in tracer units, we simply
!        use t instead of t*Hz_inv.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
!>              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>
                BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
!>
                BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)
              END DO
!
!  Impose positive definite concentrations.
!
              cff2=0.0_r8
              DO itime=1,2
                cff1=0.0_r8
                iTrcMax=idbio(1)
#ifdef IRON_LIMIT
                DO itrc=1,NBT-2
#else
                DO itrc=1,NBT
#endif
                  ibio=idbio(itrc)
                  cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                  IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                    iTrcMax=ibio
                  END IF
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
                IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                  BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
                END IF
#ifdef IRON_LIMIT
                DO itrc=NBT-1,NBT
                  ibio=idbio(itrc)
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
#endif
              END DO
!
!  Load biological tracers into local arrays.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
                Bio_old(i,k,ibio)=BioTrc(ibio,nstp)
                Bio(i,k,ibio)=BioTrc(ibio,nstp)
              END DO

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
              IF (h(i,j).le.FeHmin(ng)) THEN
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
              END IF
#endif
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
#ifdef CONST_PAR
!
!  Specify constant surface irradiance a la Powell and Spitz.
!
            PARsur(i)=158.075_r8
#else
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#endif
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Compute light attenuation as function of depth.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              IF (PARsur(i).gt.0.0_r8) THEN            ! day time
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                  Att=(AttSW(ng)+AttPhy(ng)*Bio(i,k,iPhyt))*            &
     &                (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att  ! average at cell center
                  Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
              ELSE                                     ! night time
                DO k=1,N(ng)
                  Light(i,k)=0.0_r8
                END DO
              END IF
            END DO
!
!  Phytoplankton photosynthetic growth and nitrate uptake (Vm_NO3 rate).
!  The Michaelis-Menten curve is used to describe the change in uptake
!  rate as a function of nitrate concentration. Here, PhyIS is the
!  initial slope of the P-I curve and K_NO3 is the half saturation of
!  phytoplankton nitrate uptake.
#ifdef IRON_LIMIT
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
#endif
!
            cff1=dtdays*Vm_NO3(ng)*PhyIS(ng)
            cff2=Vm_NO3(ng)*Vm_NO3(ng)
            cff3=PhyIS(ng)*PhyIS(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
#ifdef IRON_LIMIT
                FNratio=Bio(i,k,iFphy)/MAX(MinVal,Bio(i,k,iPhyt))
                FCratio=FNratio*FeN2FeC
                FCratioE=B_Fe(ng)*Bio(i,k,iFdis)**A_Fe(ng)
                Flimit=FCratio*FCratio/                                 &
     &                 (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))

                Nlimit=1.0_r8/(K_NO3(ng)+Bio(i,k,iNO3_))
                FNlim=MIN(1.0_r8,Flimit/(Bio(i,k,iNO3_)*Nlimit))
#endif
!
                cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
                cff=Bio(i,k,iPhyt)*                                     &
#ifdef IRON_LIMIT
     &              cff1*cff4*Light(i,k)*FNlim*Nlimit
#else
     &              cff1*cff4*Light(i,k)/                               &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
#endif
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff

#ifdef IRON_LIMIT
!
!  Iron uptake proportional to growth.
!
                fac=cff*Bio(i,k,iNO3_)*FNratio/                         &
     &              MAX(MinVal,Bio(i,k,iFdis))
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+fac)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*fac
!
!  Iron uptake to reach appropriate Fe:C ratio.
!
                cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
                cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
                IF (cff6.ge.0.0_r8) THEN
                  cff=cff6/MAX(MinVal,Bio(i,k,iFdis))
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)+                        &
     &                           Bio(i,k,iFdis)*cff
                ELSE
                  cff=-cff6/MAX(MinVal,Bio(i,k,iFphy))
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)+                        &
     &                           Bio(i,k,iFphy)*cff
                END IF
#endif
              END DO
            END DO
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
#ifdef IRON_LIMIT
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
#endif
!
            cff1=dtdays*ZooGR(ng)
            cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iZoop)*                                     &
     &              cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio(i,k,iPhyt)))/       &
     &              Bio(i,k,iPhyt)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/(1.0_r8+cff)
                Bio1(i,k,iZoop)=Bio(i,k,iZoop)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         Bio(i,k,iPhyt)*cff2*cff
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*ZooEEN(ng)*cff
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*ZooEED(ng)*cff
#ifdef IRON_LIMIT
                Bio1(i,k,iFphy)=Bio(i,k,iFphy)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                Bio2(i,k,iFphy)=Bio(i,k,iFphy)
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*cff*FeRR(ng)
#endif
              END DO
            END DO
!
!  Phytoplankton mortality to nutrients (PhyMRNro rate), detritus
!  (PhyMRD rate), and if applicable dissolved iron (FeRR rate).
!
            cff3=dtdays*PhyMRD(ng)
            cff2=dtdays*PhyMRN(ng)
            cff1=1.0_r8/(1.0_r8+cff2+cff3)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff1
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*cff2
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*cff3
#ifdef IRON_LIMIT
                Bio(i,k,iFphy)=Bio(i,k,iFphy)*cff1
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
#endif
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!  Zooplankton mortality to nutrients (ZooMRN rate) and Detritus
!  (ZooMRD rate).
!
              cff3=dtdays*ZooMRD(ng)
              cff2=dtdays*ZooMRN(ng)
              cff1=1.0_r8/(1.0_r8+cff2+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iZoop)*cff2
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)+                        &
     &                           Bio(i,k,iZoop)*cff3
                END DO
              END DO
!
!  Detritus breakdown to nutrients: remineralization (DetRR rate).
!
              cff2=dtdays*DetRR(ng)
              cff1=1.0_r8/(1.0_r8+cff2)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff1
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iSDet)*cff2
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Phytoplankton and Detritus
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
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
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
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
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays II.
!
!  Adjoint detritus breakdown to nutrients: remineralization
!  (DetRR rate).
!
          cff2=dtdays*DetRR(ng)
          cff1=1.0_r8/(1.0_r8+cff2)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
!>   &                          tl_Bio(i,k,iSDet)*cff2
!>
              ad_Bio(i,k,iSDet)=ad_Bio(i,k,iSDet)+                      &
     &                          cff2*ad_Bio(i,k,iNO3_)
!>            tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)*cff1
!>
              ad_Bio(i,k,iSDet)=ad_Bio(i,k,iSDet)*cff1
            END DO
          END DO
!
!  Adjoint Zooplankton mortality to nutrients (ZooMRN rate) and
!  Detritus (ZooMRD rate).
!
          cff3=dtdays*ZooMRD(ng)
          cff2=dtdays*ZooMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)+                      &
!>   &                          tl_Bio(i,k,iZoop)*cff3
!>
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+                      &
     &                          cff3*ad_Bio(i,k,iSDet)
!>            tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
!>   &                          tl_Bio(i,k,iZoop)*cff2
!>
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+                      &
     &                          cff2*ad_Bio(i,k,iNO3_)
!>            tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)*cff1
!>
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)*cff1
            END DO
          END DO
!
!  Adjoint Phytoplankton mortality to nutrients (PhyMRN rate) and
!  detritus (PhyMRD rate), and if applicable dissolved iron
!  (FeRR rate).
!
          cff3=dtdays*PhyMRD(ng)
          cff2=dtdays*PhyMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef IRON_LIMIT
!>            tl_Bio(i,k,iFdis)=tl_Bio(i,k,iFdis)+                      &
!>   &                          tl_Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
!>
              ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)+                      &
     &                          (cff2+cff3)*FeRR(ng)*ad_Bio(i,k,iFdis)
!>            tl_Bio(i,k,iFphy)=tl_Bio(i,k,iFphy)*cff1
!>
              ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)*cff1
#endif
!>            tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)+                      &
!>   &                          tl_Bio(i,k,iPhyt)*cff3
!>
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          cff3*ad_Bio(i,k,iSDet)
!>            tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
!>   &                          tl_Bio(i,k,iPhyt)*cff2
!>
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          cff2*ad_Bio(i,k,iNO3_)
!>            tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)*cff1
!>
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)*cff1
            END DO
          END DO
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
#ifdef IRON_LIMIT
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
#endif
!
          cff1=dtdays*ZooGR(ng)
          cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio1(i,k,iZoop)*                                      &
     &            cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio1(i,k,iPhyt)))/        &
     &            Bio1(i,k,iPhyt)
#ifdef IRON_LIMIT
!>            tl_Bio(i,k,iFdis)=tl_Bio(i,k,iFdis)+                      &
!>   &                          (tl_Bio(i,k,iFphy)*cff+                 &
!>   &                           Bio2(i,k,iFphy)*tl_cff)*FeRR(ng)
!>
              ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)+cff*FeRR(ng)*         &
     &                          ad_Bio(i,k,iFdis)
              ad_cff=ad_cff+Bio2(i,k,iFphy)*FeRR(ng)*ad_Bio(i,k,iFdis)
!>            tl_Bio(i,k,iFphy)=(tl_Bio(i,k,iFphy)-                     &
!>   &                           tl_cff*Bio2(i,k,iFphy))/               &
!>   &                          (1.0_r8+cff)
!>
              adfac=ad_Bio(i,k,iFphy)/(1.0_r8+cff)
              ad_cff=ad_cff-Bio2(i,k,iFphy)*adfac
              ad_Bio(i,k,iFphy)=adfac
#endif
!>            tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)+                      &
!>   &                          ZooEED(ng)*(tl_Bio(i,k,iPhyt)*cff+      &
!>   &                                      Bio(i,k,iPhyt)*tl_cff)
!>
              ad_cff=ad_cff+ZooEED(ng)*Bio(i,k,iPhyt)*ad_Bio(i,k,iSDet)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          ZooEED(ng)*cff*ad_Bio(i,k,iSDet)
!>            tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
!>   &                          ZooEEN(ng)*(tl_Bio(i,k,iPhyt)*cff+      &
!>   &                                      Bio(i,k,iPhyt)*tl_cff)
!>
              ad_cff=ad_cff+                                            &
     &               ZooEEN(ng)*Bio(i,k,iPhyt)*ad_Bio(i,k,iNO3_)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          ZooEEN(ng)*cff*ad_Bio(i,k,iNO3_)
!>            tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)+                      &
!>   &                          cff2*(tl_Bio(i,k,iPhyt)*cff+            &
!>   &                                Bio(i,k,iPhyt)*tl_cff)
!>
              ad_cff=ad_cff+                                            &
     &               cff2*Bio(i,k,iPhyt)*ad_Bio(i,k,iZoop)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          cff2*cff*ad_Bio(i,k,iZoop)
!>            tl_Bio(i,k,iPhyt)=(tl_Bio(i,k,iPhyt)-                     &
!>   &                           tl_cff*Bio(i,k,iPhyt))/                &
!>   &                          (1.0_r8+cff)
!>
              adfac=ad_Bio(i,k,iPhyt)/(1.0_r8+cff)
              ad_cff=ad_cff-Bio(i,k,iPhyt)*adfac
              ad_Bio(i,k,iPhyt)=adfac
!>            tl_cff=(tl_Bio(i,k,iZoop)*                                &
!>   &                cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio1(i,k,iPhyt)))+    &
!>   &                Bio1(i,k,iZoop)*Ivlev(ng)*tl_Bio(i,k,iPhyt)*cff1* &
!>   &                EXP(-Ivlev(ng)*Bio1(i,k,iPhyt))-                  &
!>   &                tl_Bio(i,k,iPhyt)*cff)/                           &
!>   &               Bio1(i,k,iPhyt)
!>
              fac=EXP(-Ivlev(ng)*Bio1(i,k,iPhyt))
              adfac=ad_cff/Bio1(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)-adfac*cff
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          Bio1(i,k,iZoop)*Ivlev(ng)*              &
     &                          adfac*cff1*fac
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+                      &
     &                          adfac*cff1*(1.0_r8-fac)
              ad_cff=0.0_r8
            END DO
          END DO
!
!  Compute appropriate basic state arrays I.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
!  NOTE: In the following code, t(:,:,:,nnew,:) should be in units of
!        tracer times depth. However the basic state (nstp and nnew
!        indices) that is read from the forward file is in units of
!        tracer. Since BioTrc(ibio,nnew) is in tracer units, we simply
!        use t instead of t*Hz_inv.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
!>              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>
                BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
!>              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
!>
                BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)
              END DO
!
!  Impose positive definite concentrations.
!
              cff2=0.0_r8
              DO itime=1,2
                cff1=0.0_r8
                iTrcMax=idbio(1)
#ifdef IRON_LIMIT
                DO itrc=1,NBT-2
#else
                DO itrc=1,NBT
#endif
                  ibio=idbio(itrc)
                  cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                  IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                    iTrcMax=ibio
                  END IF
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
                IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                  BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
                END IF
#ifdef IRON_LIMIT
                DO itrc=NBT-1,NBT
                  ibio=idbio(itrc)
                  BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
                END DO
#endif
              END DO
!
!  Load biological tracers into local arrays.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
                Bio_old(i,k,ibio)=BioTrc(ibio,nstp)
                Bio(i,k,ibio)=BioTrc(ibio,nstp)
              END DO

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
              IF (h(i,j).le.FeHmin(ng)) THEN
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
              END IF
#endif
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
#ifdef CONST_PAR
!
!  Specify constant surface irradiance a la Powell and Spitz.
!
            PARsur(i)=158.075_r8
#else
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#endif
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Compute light attenuation as function of depth.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              IF (PARsur(i).gt.0.0_r8) THEN            ! day time
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                  Att=(AttSW(ng)+AttPhy(ng)*Bio(i,k,iPhyt))*            &
     &                (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att  ! average at cell center
                  Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
              ELSE                                     ! night time
                DO k=1,N(ng)
                  Light(i,k)=0.0_r8
                END DO
              END IF
            END DO
!
!  Phytoplankton photosynthetic growth and nitrate uptake (Vm_NO3 rate).
!  The Michaelis-Menten curve is used to describe the change in uptake
!  rate as a function of nitrate concentration. Here, PhyIS is the
!  initial slope of the P-I curve and K_NO3 is the half saturation of
!  phytoplankton nitrate uptake.
#ifdef IRON_LIMIT
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
#endif
!
            cff1=dtdays*Vm_NO3(ng)*PhyIS(ng)
            cff2=Vm_NO3(ng)*Vm_NO3(ng)
            cff3=PhyIS(ng)*PhyIS(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
#ifdef IRON_LIMIT
!
!  Calculate growth reduction factor due to iron limitation.
!
                FNratio=Bio(i,k,iFphy)/MAX(MinVal,Bio(i,k,iPhyt))
                FCratio=FNratio*FeN2FeC
                FCratioE=B_Fe(ng)*Bio(i,k,iFdis)**A_Fe(ng)
                Flimit=FCratio*FCratio/                                 &
     &                 (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))

                Nlimit=1.0_r8/(K_NO3(ng)+Bio(i,k,iNO3_))
                FNlim=MIN(1.0_r8,Flimit/(Bio(i,k,iNO3_)*Nlimit))
#endif
                cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
                cff=Bio(i,k,iPhyt)*                                     &
#ifdef IRON_LIMIT
     &              cff1*cff4*Light(i,k)*FNlim*Nlimit
#else
     &              cff1*cff4*Light(i,k)/                               &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
#endif
                Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff

#ifdef IRON_LIMIT
!
!  Iron uptake proportional to growth.
!
                fac=cff*Bio(i,k,iNO3_)*FNratio/                         &
     &              MAX(MinVal,Bio(i,k,iFdis))
                Bio1(i,k,iFdis)=Bio(i,k,iFdis)
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+fac)
                Bio2(i,k,iFdis)=Bio(i,k,iFdis)
                Bio1(i,k,iFphy)=Bio(i,k,iFphy)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*fac
                Bio2(i,k,iFphy)=Bio(i,k,iFphy)
!
!  Iron uptake to reach appropriate Fe:C ratio.
!
                cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
                cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
                IF (cff6.ge.0.0_r8) then
                  cff=cff6/MAX(MinVal,Bio(i,k,iFdis))
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)+                        &
     &                           Bio(i,k,iFdis)*cff
                ELSE
                  cff=-cff6/MAX(MinVal,Bio(i,k,iFphy))
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)+                        &
     &                           Bio(i,k,iFphy)*cff
                END IF
#endif
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
#ifdef IRON_LIMIT
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
#endif
!
              cff1=dtdays*ZooGR(ng)
              cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=Bio(i,k,iZoop)*                                   &
     &                cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio(i,k,iPhyt)))/     &
     &                Bio(i,k,iPhyt)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/(1.0_r8+cff)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)+                        &
     &                           Bio(i,k,iPhyt)*cff2*cff
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iPhyt)*ZooEEN(ng)*cff
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)+                        &
     &                           Bio(i,k,iPhyt)*ZooEED(ng)*cff
#ifdef IRON_LIMIT
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)+                        &
     &                           Bio(i,k,iFphy)*cff*FeRR(ng)
#endif
                END DO
              END DO
!
!  Phytoplankton mortality to nutrients (PhyMRNro rate), detritus
!  (PhyMRD rate), and if applicable dissolved iron (FeRR rate).
!
              cff3=dtdays*PhyMRD(ng)
              cff2=dtdays*PhyMRN(ng)
              cff1=1.0_r8/(1.0_r8+cff2+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff1
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iPhyt)*cff2
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)+                        &
     &                           Bio(i,k,iPhyt)*cff3
#ifdef IRON_LIMIT
                  Bio(i,k,iFphy)=Bio(i,k,iFphy)*cff1
                  Bio(i,k,iFdis)=Bio(i,k,iFdis)+                        &
     &                           Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
#endif
                END DO
              END DO
!
!  Zooplankton mortality to nutrients (ZooMRN rate) and Detritus
!  (ZooMRD rate).
!
              cff3=dtdays*ZooMRD(ng)
              cff2=dtdays*ZooMRN(ng)
              cff1=1.0_r8/(1.0_r8+cff2+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iZoop)*cff2
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)+                        &
     &                           Bio(i,k,iZoop)*cff3
                END DO
              END DO
!
!  Detritus breakdown to nutrients: remineralization (DetRR rate).
!
              cff2=dtdays*DetRR(ng)
              cff1=1.0_r8/(1.0_r8+cff2)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff1
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iSDet)*cff2
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Phytoplankton and Detritus
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
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
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
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
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays I.
!
!  Adjoint Phytoplankton photosynthetic growth and nitrate uptake
!  (Vm_NO3 rate). The Michaelis-Menten curve is used to describe the
!  change in uptake rate as a function of nitrate concentration.
!  Here, PhyIS is the initial slope of the P-I curve and K_NO3 is the
!  half saturation of phytoplankton nitrate uptake.
#ifdef IRON_LIMIT
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
#endif
!
          cff1=dtdays*Vm_NO3(ng)*PhyIS(ng)
          cff2=Vm_NO3(ng)*Vm_NO3(ng)
          cff3=PhyIS(ng)*PhyIS(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef IRON_LIMIT
!
!  Adjoint of iron uptake to reach appropriate Fe:C ratio.
!
              FNratio=Bio1(i,k,iFphy)/MAX(MinVal,Bio1(i,k,iPhyt))
              FCratio=FNratio*FeN2FeC
              FCratioE=B_Fe(ng)*Bio1(i,k,iFdis)**A_Fe(ng)
              Flimit=FCratio*FCratio/                                   &
     &               (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))

              Nlimit=1.0_r8/(K_NO3(ng)+Bio1(i,k,iNO3_))
              FNlim=MIN(1.0_r8,Flimit/(Bio1(i,k,iNO3_)*Nlimit))

              cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
              cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
              IF (cff6.ge.0.0_r8) THEN
                fac1=MAX(MinVal,Bio2(i,k,iFdis))
                cff=cff6/fac1
!>              tl_Bio(i,k,iFphy)=tl_Bio(i,k,iFphy)+                    &
!>   &                            tl_Bio(i,k,iFdis)*cff+                &
!>   &                            Bio(i,k,iFdis)*tl_cff
!>
                ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)+                    &
     &                            cff*ad_Bio(i,k,iFphy)
                ad_cff=ad_cff+Bio(i,k,iFdis)*ad_Bio(i,k,iFphy)
!>              tl_Bio(i,k,iFdis)=(tl_Bio(i,k,iFdis)-                   &
!>   &                             tl_cff*Bio(i,k,iFdis))/              &
!>   &                            (1.0_r8+cff)
!>
                adfac=ad_Bio(i,k,iFdis)/(1.0_r8+cff)
                ad_cff=ad_cff-Bio(i,k,iFdis)*adfac
                ad_Bio(i,k,iFdis)=adfac
!>              tl_cff=(tl_cff6-tl_fac1*cff)/fac1
!>
                adfac=ad_cff/fac1
                ad_cff6=ad_cff6+adfac
                ad_fac1=ad_fac1-cff*adfac
                ad_cff=0.0_r8
!>              tl_fac1=(0.5_r8-SIGN(0.5_r8,MinVal-Bio2(i,k,iFdis)))*   &
!>   &                  tl_Bio(i,k,iFdis)
!>
                ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)+                    &
     &                            (0.5_r8-                              &
     &                             SIGN(0.5_r8,                         &
     &                                  MinVal-Bio2(i,k,iFdis)))*ad_fac1
                ad_fac1=0.0_r8
              ELSE
                fac1=-MAX(MinVal,Bio2(i,k,iFphy))
                cff=cff6/fac1
!>              tl_Bio(i,k,iFdis)=tl_Bio(i,k,iFdis)+                    &
!>   &                            tl_Bio(i,k,iFphy)*cff+                &
!>   &                            Bio(i,k,iFphy)*tl_cff
!>
                ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)+                    &
     &                            cff*ad_Bio(i,k,iFdis)
                ad_cff=ad_cff+Bio(i,k,iFphy)*ad_Bio(i,k,iFdis)
!>              tl_Bio(i,k,iFphy)=(tl_Bio(i,k,iFphy)-                   &
!>   &                             tl_cff*Bio(i,k,iFphy))/              &
!>   &                            (1.0_r8+cff)
!>
                adfac=ad_Bio(i,k,iFphy)/(1.0_r8+cff)
                ad_cff=ad_cff-Bio(i,k,iFphy)*adfac
                ad_Bio(i,k,iFphy)=adfac
!>              tl_cff=(tl_cff6-tl_fac1*cff)/fac1
!>
                adfac=ad_cff/fac1
                ad_cff6=ad_cff6+adfac
                ad_fac1=ad_fac1-cff*adfac
                ad_cff=0.0_r8
!>              tl_fac1=-(0.5_r8-SIGN(0.5_r8,MinVal-Bio2(i,k,iFphy)))*  &
!>   &                  tl_Bio(i,k,iFphy)
!>
                ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)-                    &
     &                            (0.5_r8-                              &
     &                             SIGN(0.5_r8,                         &
     &                                  MinVal-Bio2(i,k,iFphy)))*ad_fac1
                ad_fac1=0.0_r8
              END IF
!>            tl_cff6=(tl_Bio(i,k,iPhyt)*cff5+                          &
!>   &                 Bio(i,k,iPhyt)*tl_cff5)*FeC2FeN
!>
              adfac=ad_cff6*FeC2FeN
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+cff5*adfac
              ad_cff5=ad_cff5+Bio(i,k,iPhyt)*adfac
              ad_cff6=0.0_r8
!>            tl_cff5=dtdays*(tl_FCratioE-tl_FCratio)/T_Fe(ng)
!>
              adfac=dtdays*ad_cff5/T_Fe(ng)
              ad_FCratioE=ad_FCratioE+adfac
              ad_FCratio=ad_FCratio-adfac
              ad_cff5=0.0_r8
#endif
              cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
              cff=Bio1(i,k,iPhyt)*                                      &
#ifdef IRON_LIMIT
     &            cff1*cff4*Light(i,k)*FNlim*Nlimit
#else
     &            cff1*cff4*Light(i,k)/                                 &
     &            (K_NO3(ng)+Bio1(i,k,iNO3_))
#endif
#ifdef IRON_LIMIT
!
!  Iron uptake proportional to growth.
!
              fac1=MAX(MinVal,Bio1(i,k,iFdis))
              fac2=1.0_r8/fac1
              fac=cff*Bio(i,k,iNO3_)*FNratio*fac2
!>            tl_Bio(i,k,iFphy)=tl_Bio(i,k,iFphy)+                      &
!>   &                          tl_Bio(i,k,iFdis)*fac+                  &
!>   &                          Bio2(i,k,iFdis)*tl_fac
!>
              ad_fac=ad_fac+Bio2(i,k,iFdis)*ad_Bio(i,k,iFphy)
              ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)+                      &
     &                          fac*ad_Bio(i,k,iFphy)
!>            tl_Bio(i,k,iFdis)=(tl_Bio(i,k,iFdis)-                     &
!>   &                           tl_fac*Bio2(i,k,iFdis))/               &
!>   &                          (1.0_r8+fac)
!>
              adfac=ad_Bio(i,k,iFdis)/(1.0_r8+fac)
              ad_fac=ad_fac-Bio2(i,k,iFdis)*adfac
              ad_Bio(i,k,iFdis)=adfac
!>            tl_fac=FNratio*fac2*(tl_cff*Bio(i,k,iNO3_)+               &
!>   &                             cff*ad_Bio(i,k,iNO3_))+              &
!>   &               cff*Bio(i,k,iNO3_)*(tl_FNratio*fac2+               &
!>   &                                   FNratio*tl_fac2)
!>
              adfac1=FNratio*fac2*ad_fac
              adfac2=cff*Bio(i,k,iNO3_)*ad_fac
              ad_cff=ad_cff+Bio(i,k,iNO3_)*adfac1
              ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+cff*adfac1
              ad_FNratio=ad_FNratio+fac2*adfac2
              ad_fac2=ad_fac2+FNratio*adfac2
              ad_fac=0.0_r8
!>            tl_fac2=-fac2*fac2*tl_fac1
!>
              ad_fac1=ad_fac1-fac2*fac2*ad_fac2
              ad_fac2=0.0_r8
!>            tl_fac1=(0.5_r8-SIGN(0.5_r8,MinVal-Bio1(i,k,iFdis)))*     &
!>   &                tl_Bio(i,k,iFdis)
!>
              ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)+                      &
     &                          (0.5_r8-                                &
     &                           SIGN(0.5_r8,MinVal-Bio1(i,k,iFdis)))*  &
     &                          ad_fac1
              ad_fac1=0.0_r8
#endif
!
!  Adjoint of phytoplankton photosynthetic growth and nitrate uptake.
!
!>            tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+                      &
!>   &                          tl_Bio(i,k,iNO3_)*cff+                  &
!>   &                          Bio(i,k,iNO3_)*tl_cff
!>
              ad_cff=ad_cff+Bio(i,k,iNO3_)*ad_Bio(i,k,iPhyt)
              ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+                      &
     &                          cff*ad_Bio(i,k,iPhyt)
!>            tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-                     &
!>   &                           tl_cff*Bio(i,k,iNO3_))/                &
!>   &                          (1.0_r8+cff)
!>
              adfac=ad_Bio(i,k,iNO3_)/(1.0_r8+cff)
              ad_cff=ad_cff-Bio(i,k,iNO3_)*adfac
              ad_Bio(i,k,iNO3_)=adfac

#ifdef IRON_LIMIT
!>            tl_cff=tl_Bio(i,k,iPhyt)*                                 &
!>   &               cff1*cff4*Light(i,k)*FNlim*Nlimit+                 &
!>   &               Bio1(i,k,iPhyt)*cff1*cff4*                         &
!>   &               (tl_Light(i,k)*FNlim*Nlimit+                       &
!>   &                Light(i,k)*tl_FNlim*Nlimit+                       &
!>   &                Light(i,k)*FNlim*tl_Nlimit)+                      &
!>   &               Bio1(i,k,iPhyt)*cff1*tl_cff4*                      &
!>   &               Light(i,k)*FNlim*Nlimit
!>
              adfac1=cff1*cff4*ad_cff
              adfac2=adfac1*Bio1(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          Light(i,k)*FNlim*Nlimit*adfac1
              ad_Light(i,k)=ad_Light(i,k)+                              &
     &                      FNlim*Nlimit*adfac2
              ad_FNlim=ad_FNlim+                                        &
     &                 Light(i,k)*Nlimit*adfac2
              ad_Nlimit=ad_Nlimit+                                      &
     &                  Light(i,k)*FNlim*adfac2
              ad_cff4=ad_cff4+                                          &
     &                Bio1(i,k,iPhyt)*cff1*Light(i,k)*FNlim*Nlimit*     &
     &                ad_cff
              ad_cff=0.0_r8
#else
!>            tl_cff=(tl_Bio(i,k,iPhyt)*cff1*cff4*Light(i,k)+           &
!>   &                Bio1(i,k,iPhyt)*cff1*                             &
!>   &                (tl_cff4*Light(i,k)+cff4*tl_Light(i,k))-          &
!>   &                tl_Bio(i,k,iNO3_)*cff)/                           &
!>   &               (K_NO3(ng)+Bio1(i,k,iNO3_))
!>
              adfac=ad_cff/(K_NO3(ng)+Bio1(i,k,iNO3_))
              adfac1=adfac*Bio1(i,k,iPhyt)*cff1
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          cff1*cff4*Light(i,k)*adfac
              ad_cff4=adfac1*Light(i,k)
              ad_Light(i,k)=ad_Light(i,k)+                              &
     &                      adfac1*cff4
              ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)-                      &
     &                          adfac*cff
              ad_cff=0.0_r8
#endif
!>            tl_cff4=-cff3*tl_Light(i,k)*Light(i,k)*cff4*cff4*cff4
!>
              ad_Light(i,k)=ad_Light(i,k)-                              &
     &                      cff3*Light(i,k)*                            &
     &                      cff4*cff4*cff4*ad_cff4
              ad_cff4=0.0_r8

#ifdef IRON_LIMIT
!
!  Adjoint of  calculate growth reduction factor due to iron limitation.
!
              Nlimit=1.0_r8/(K_NO3(ng)+Bio1(i,k,iNO3_))
              fac1=Flimit/(Bio1(i,k,iNO3_)*Nlimit)
!>            tl_FNlim=(0.5_r8+SIGN(0.5_r8,1.0_r8-fac1))*tl_fac1
!>
              ad_fac1=(0.5_r8+SIGN(0.5_r8,1.0_r8-fac1))*ad_FNlim
              ad_FNlim=0.0_r8
!>            tl_fac1=tl_Flimit/(Bio1(i,k,iNO3_)*Nlimit)-               &
!>   &                (tl_Bio(i,k,iNO3_)*Nlimit+                        &
!>   &                 Bio1(i,k,iNO3_)*tl_Nlimit)*fac1/                 &
!>   &                (Bio1(i,k,iNO3_)*Nlimit)
!>
              adfac1=ad_fac1/(Bio1(i,k,iNO3_)*Nlimit)
              adfac2=adfac1*fac1
              ad_Flimit=ad_Flimit+adfac1
              ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)-Nlimit*adfac2
              ad_Nlimit=ad_Nlimit-Bio1(i,k,iNO3_)*adfac2
              ad_fac1=0.0_r8
!>            tl_Nlimit=-tl_Bio(i,k,iNO3_)*Nlimit*Nlimit
!>
              ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)-                      &
     &                          ad_Nlimit*Nlimit*Nlimit
              ad_Nlimit=0.0_r8
!>            tl_Flimit=2.0_r8*(tl_FCratio*FCratio-                     &
!>   &                          tl_FCratio*FCratio*Flimit)/             &
!>   &                  (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))
!>
              adfac=2.0_r8*ad_Flimit/                                   &
     &              (FCratio*FCratio+K_FeC(ng)*K_FeC(ng))
              ad_FCratio=ad_FCratio+                                    &
     &                   FCratio*adfac-                                 &
     &                   FCratio*Flimit*adfac
              ad_Flimit=0.0_r8
!>            tl_FCratioE=A_Fe(ng)*B_Fe(ng)*                            &
!>   &                    Bio1(i,k,iFdis)**(A_Fe(ng)-1.0_r8)*           &
!>   &                    tl_Bio(i,k,iFdis)
!>
              ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)+                      &
     &                          A_Fe(ng)*B_Fe(ng)*                      &
     &                          Bio1(i,k,iFdis)**(A_Fe(ng)-1.0_r8)*     &
     &                          ad_FCratioE
              ad_FCratioE=0.0_r8
!>            tl_FCratio=tl_FNratio*FeN2FeC
!>
              ad_FNratio=ad_FNratio+FeN2FeC*ad_FCratio
              ad_FCratio=0.0_r8

              fac1=MAX(MinVal,Bio1(i,k,iPhyt))
!>            tl_FNratio=(tl_Bio(i,k,iFphy)-tl_fac1*FNratio)/fac1
!>
              adfac=ad_FNratio/fac1
              ad_Bio(i,k,iFphy)=ad_Bio(i,k,iFphy)+adfac
              ad_fac1=ad_fac1-ad_FNratio*FNratio/fac1
              ad_FNratio=0.0_r8
!>            tl_fac1=(0.5_r8-SIGN(0.5_r8,MinVal-Bio1(i,k,iPhyt)))*     &
!>   &                tl_Bio(i,k,iPhyt)
!>
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          (0.5_r8-                                &
     &                           SIGN(0.5_r8,MinVal-Bio1(i,k,iPhyt)))*  &
     &                          ad_fac1
              ad_fac1=0.0_r8
#endif
            END DO
          END DO
!
!  Compute adjoint light attenuation as function of depth.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=1,N(ng)
!
! Compute the basic state PAR appropriate for each level.
!
                PAR=PARsur(i)
                DO kk=N(ng),k,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                  Att=(AttSW(ng)+AttPhy(ng)*Bio1(i,kk,iPhyt))*          &
     &                (z_w(i,j,kk)-z_w(i,j,kk-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att  ! average at cell center
                  PAR1=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!
!  Adjoint of light attenuation at the bottom of the grid cell. It is
!  the starting PAR value for the next (deeper) vertical grid cell.
!
!>              tl_PAR=tl_Itop*ExpAtt+Itop*tl_ExpAtt
!>
                ad_ExpAtt=ad_ExpAtt+Itop*ad_PAR
                ad_Itop=ad_Itop+ExpAtt*ad_PAR
                ad_PAR=0.0_r8
!
!  Adjoint of compute average light attenuation for each grid cell.
!  Here, AttSW is the light attenuation due to seawater and AttPhy is
!  the attenuation due to phytoplankton (self-shading coefficient).
!
!>              tl_Light(i,k)=tl_PAR
!>
                ad_PAR=ad_PAR+ad_Light(i,k)
                ad_Light(i,k)=0.0_r8
!>              tl_PAR=(-tl_Att*PAR1+tl_Itop*(1.0_r8-ExpAtt)-           &
!>   &                  Itop*tl_ExpAtt)/Att
!>
                adfac=ad_PAR/Att
                ad_Att=ad_Att-PAR1*adfac
                ad_ExpAtt=ad_ExpAtt-Itop*adfac
                ad_Itop=ad_Itop+(1.0_r8-ExpAtt)*adfac
                ad_PAR=0.0_r8
!>              tl_Itop=tl_PAR
!>
                ad_PAR=ad_PAR+ad_Itop
                ad_Itop=0.0_r8
!>              tl_ExpAtt=-ExpAtt*tl_Att
!>
                ad_Att=ad_Att-ExpAtt*ad_ExpAtt
                ad_ExpAtt=0.0_r8
!>              tl_Att=AttPhy(ng)*tl_Bio(i,k,iPhyt)*                    &
!>   &                 (z_w(i,j,k)-z_w(i,j,k-1))+                       &
!>   &                 (AttSW(ng)+AttPhy(ng)*Bio1(i,k,iPhyt))*          &
!>   &                 (tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
!>
                adfac=(AttSW(ng)+AttPhy(ng)*Bio1(i,k,iPhyt))*ad_Att
                ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                    &
     &                            AttPhy(ng)*(z_w(i,j,k)-z_w(i,j,k-1))* &
     &                            ad_Att
                ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-adfac
                ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+adfac
                ad_Att=0.0_r8
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
!>              tl_Light(i,k)=0.0_r8
!>
                ad_Light(i,k)=0.0_r8
              END DO
            END IF
!>          tl_PAR=tl_PARsur(i)
!>
            ad_PARsur(i)=ad_PARsur(i)+ad_PAR
            ad_PAR=0.0_r8
          END DO

        END DO ITER_LOOP1
!
!  Calculate adjoint surface Photosynthetically Available Radiation
!  (PAR).  The net shortwave radiation is scaled back to Watts/m2
!  and multiplied by the fraction that is photosynthetically
!  available, PARfrac.
!
        DO i=Istr,Iend
#ifdef CONST_PAR
!
!  Specify constant surface irradiance a la Powell and Spitz.
!
!>        tl_PARsur(i)=0.0_r8
!>
!!        ad_PARsur(i)=0.0_r8
#else
!>        tl_PARsur(i)=(tl_PARfrac(ng)*srflx(i,j)+                      &
!>   &                  PARfrac(ng)*tl_srflx(i,j))*rho0*Cp
!>
          adfac=rho0*Cp*ad_PARsur(i)
          ad_srflx(i,j)=ad_srflx(i,j)+PARfrac(ng)*adfac
          ad_PARfrac(ng)=ad_PARfrac(ng)+srflx(i,j)*adfac
!!        ad_PARsur(i)=0.0_r8
#endif
        END DO
!
!  Restrict biological tracer to be positive definite. If a negative
!  concentration is detected, nitrogen is drawn from the most abundant
!  pool to supplement the negative pools to a lower limit of MinVal
!  which is set to 1E-6 above.
!
        DO k=1,N(ng)
          DO i=Istr,Iend

#if defined IRON_LIMIT && defined IRON_RELAX
!
!  Adjoint of relax dissolved iron at coast (h <= FeHim) to a constant
!  value (FeMax) over a time scale (FeNudgTime; days) to simulate
!  sources at the shelf.
!
            IF (h(i,j).le.FeHmin(ng)) THEN
!>            tl_Bio(i,k,iFdis)=tl_Bio(i,k,iFdis)-                      &
!>   &                          FeNudgCoef*tl_Bio(i,k,iFdis)
!>
              ad_Bio(i,k,iFdis)=ad_Bio(i,k,iFdis)-                      &
     &                          FeNudgCoef*ad_Bio(i,k,iFdis)
            END IF
#endif
!
!  Adjoint load biological tracers into local arrays.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
!>            tl_Bio(i,k,ibio)=tl_BioTrc(ibio,nstp)
!>
              ad_BioTrc(ibio,nstp)=ad_BioTrc(ibio,nstp)+                &
     &                             ad_Bio(i,k,ibio)
              ad_Bio(i,k,ibio)=0.0_r8
!>            tl_Bio_old(i,k,ibio)=tl_BioTrc(ibio,nstp)
!>
              ad_BioTrc(ibio,nstp)=ad_BioTrc(ibio,nstp)+                &
     &                             ad_Bio_old(i,k,ibio)
              ad_Bio_old(i,k,ibio)=0.0_r8
            END DO
!
!  Adjoint positive definite concentrations.
!
            DO itime=1,2
              DO itrc=1,NBT
                ibio=idbio(itrc)
!
!  The basic state (nstp and nnew indices) that is read from the
!  forward file is in units of tracer. Since BioTrc(ibio,:) is in
!  tracer units, we simply use t instead of t*Hz_inv.
!
                BioTrc(ibio,itime)=t(i,j,k,itime,ibio)
                BioTrc1(ibio,itime)=BioTrc(ibio,itime)
              END DO
            END DO
!
            cff2=0.0_r8
            DO itime=1,2
              cff1=0.0_r8
              iTrcMax=idbio(1)
#ifdef IRON_LIMIT
              DO itrc=1,NBT-2
#else
              DO itrc=1,NBT
#endif
                ibio=idbio(itrc)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                  iTrcMax=ibio
                END IF
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
              IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
              END IF

#ifdef IRON_LIMIT
              DO itrc=NBT-1,NBT
                ibio=idbio(itrc)
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc1(ibio,itime))
!>              tl_BioTrc(ibio,itime)=(0.5_r8-                          &
!>   &                                 SIGN(0.5_r8,                     &
!>   &                                      MinVal-                     &
!>   &                                      BioTrc1(ibio,itime)))*      &
!>   &                                tl_BioTrc(ibio,itime)
!>
                ad_BioTrc(ibio,itime)=(0.5_r8-                          &
     &                                 SIGN(0.5_r8,                     &
     &                                      MinVal-                     &
     &                                      BioTrc1(ibio,itime)))*      &
     &                                ad_BioTrc(ibio,itime)
              END DO
#endif
!
!  Recompute BioTrc again because iTrcMax has changed.
!
              cff1=0.0_r8
              iTrcMax=idbio(1)
#ifdef IRON_LIMIT
              DO itrc=1,NBT-2
#else
              DO itrc=1,NBT
#endif
                ibio=idbio(itrc)
!
!  The basic state (nstp and nnew indices) that is read from the
!  forward file is in units of tracer. Since BioTrc(ibio,:) is in
!  tracer units, we simply use t instead of t*Hz_inv.
!
                BioTrc(ibio,itime)=t(i,j,k,itime,ibio)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                  iTrcMax=ibio
                END IF
                BioTrc1(ibio,itime)=BioTrc(ibio,itime)
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
              IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
!>              tl_BioTrc(iTrcMax,itime)=tl_BioTrc(iTrcMax,itime)-      &
!>   &                                   tl_cff1
!>
                ad_cff1=-ad_BioTrc(iTrcMax,itime)
              END IF

              cff1=0.0_r8
#ifdef IRON_LIMIT
              DO itrc=1,NBT-2
#else
              DO itrc=1,NBT
#endif
                ibio=idbio(itrc)
!
!  The basic state (nstp and nnew indices) that is read from the
!  forward file is in units of tracer. Since BioTrc(ibio,:) is in
!  tracer units, we simply use t instead of t*Hz_inv.
!
                BioTrc(ibio,itime)=t(i,j,k,itime,ibio)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
!>              tl_BioTrc(ibio,itime)=(0.5_r8-                          &
!>   &                                 SIGN(0.5_r8,                     &
!>   &                                      MinVal-                     &
!>   &                                      BioTrc1(ibio,itime)))*      &
!>   &                                tl_BioTrc(ibio,itime)
!>
                ad_BioTrc(ibio,itime)=(0.5_r8-                          &
     &                                 SIGN(0.5_r8,                     &
     &                                      MinVal-                     &
     &                                      BioTrc1(ibio,itime)))*      &
     &                                ad_BioTrc(ibio,itime)
!>              tl_cff1=tl_cff1-                                        &
!>   &                  (0.5_r8-SIGN(0.5_r8,                            &
!>   &                               BioTrc(ibio,itime)-MinVal))*       &
!>   &                  tl_BioTrc(ibio,itime)
!>
                ad_BioTrc(ibio,itime)=ad_BioTrc(ibio,itime)-            &
     &                                (0.5_r8-SIGN(0.5_r8,              &
     &                                             BioTrc(ibio,itime)-  &
     &                                             MinVal))*ad_cff1
              END DO
              ad_cff1=0.0_r8
            END DO
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
!  NOTE: In the following code, t(:,:,:,nnew,:) should be in units of
!        tracer times depth. However the basic state (nstp and nnew
!        indices) that is read from the forward file is in units of
!        tracer. Since BioTrc(ibio,nnew) is in tracer units, we simply
!        use t instead of t*Hz_inv.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
!>            tl_BioTrc(ibio,nnew)=tl_t(i,j,k,nnew,ibio)*               &
!>   &                             Hz_inv(i,k)+                         &
!>   &                             t(i,j,k,nnew,ibio)*Hz(i,j,k)*        &
!>   &                             tl_Hz_inv(i,k)
!>
              ad_Hz_inv(i,k)=ad_Hz_inv(i,k)+                            &
     &                       t(i,j,k,nnew,ibio)*Hz(i,j,k)*              &
     &                       ad_BioTrc(ibio,nnew)
              ad_t(i,j,k,nnew,ibio)=ad_t(i,j,k,nnew,ibio)+              &
     &                              Hz_inv(i,k)*ad_BioTrc(ibio,nnew)
              ad_BioTrc(ibio,nnew)=0.0_r8
!>            tl_BioTrc(ibio,nstp)=tl_t(i,j,k,nstp,ibio)
!>
              ad_t(i,j,k,nstp,ibio)=ad_t(i,j,k,nstp,ibio)+              &
     &                              ad_BioTrc(ibio,nstp)
              ad_BioTrc(ibio,nstp)=0.0_r8
            END DO
          END DO
        END DO
!
!  Adjoint inverse thickness to avoid repeated divisions.
!
        DO k=2,N(ng)-1
          DO i=Istr,Iend
!>          tl_Hz_inv3(i,k)=-Hz_inv3(i,k)*Hz_inv3(i,k)*                 &
!>   &                      (tl_Hz(i,j,k-1)+tl_Hz(i,j,k)+               &
!>   &                       tl_Hz(i,j,k+1))
!>
            adfac=Hz_inv3(i,k)*Hz_inv3(i,k)*ad_Hz_inv3(i,k)
            ad_Hz(i,j,k-1)=ad_Hz(i,j,k-1)-adfac
            ad_Hz(i,j,k  )=ad_Hz(i,j,k  )-adfac
            ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)-adfac
            ad_Hz_inv3(i,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!>          tl_Hz_inv2(i,k)=-Hz_inv2(i,k)*Hz_inv2(i,k)*                 &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i,j,k+1))
!>
            adfac=Hz_inv2(i,k)*Hz_inv2(i,k)*ad_Hz_inv2(i,k)
            ad_Hz(i,j,k  )=ad_Hz(i,j,k  )-adfac
            ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)-adfac
            ad_Hz_inv2(i,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)
!>
            ad_Hz(i,j,k)=ad_Hz(i,j,k)-                                  &
     &                   Hz_inv(i,k)*Hz_inv(i,k)*ad_Hz_inv(i,k)
            ad_Hz_inv(i,k)=0.0_r8
          END DO
        END DO

      END DO J_LOOP
!
!  Set adjoint vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
!>    tl_Wbio(2)=tl_wDet(ng)          ! Small detritus
!>
      ad_wDet(ng)=ad_wDet(ng)+ad_Wbio(2)
      ad_Wbio(2)=0.0_r8
!>    tl_Wbio(1)=tl_wPhy(ng)          ! Phytoplankton
!>
      ad_wPhy(ng)=ad_wPhy(ng)+ad_Wbio(1)
      ad_Wbio(1)=0.0_r8

      RETURN
      END SUBROUTINE ad_biology_tile
