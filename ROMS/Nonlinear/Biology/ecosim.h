      SUBROUTINE biology (ng, tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!*************************************************** W. Paul Bissett ***
!  Copyright (c) 1997 W. Paul Bissett, FERI                            !
!***********************************************************************
!                                                                      !
!  The EcoSim code has been developed for research purposes only. It   !
!  consists of unpublished, proprietary formulations protected under   !
!  U.S. copyright law. It is freely available on request from the      !
!  Florida Environmental Research Institute (FERI). Commercial usage   !
!  of these formulations is forbidden without express written          !
!  permission from FERI. All rights reserved.                          !
!                                                                      !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This routine computes the EcoSim sources and sinks and adds them    !
!  to the global biological fields.                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Bissett, W.P., J.J. Walsh, D.A. Dieterle, K.L. Carder, 1999:      !
!      Carbon cycling in the upper waters of the Sargasso Sea: I.      !
!      Numerical  simulation of  differential carbon and nitrogen      !
!      fluxes,  Deep-Sea Res., 46, 205-269.                            !
!                                                                      !
!    Bissett, W.P., K.L. Carder, J.J. Walsh, D.A. Dieterle, 1999:      !
!      Carbon cycling in the upper waters of the Sargasso Sea: II.     !
!      Numerical  simulation  of  apparent  and  inherent optical      !
!      properties, Deep-Sea Res., 46, 271-317                          !
!                                                                      !
!  NOTES to EcoSim:                                                    !
!                                                                      !
!  * This version uses a descending index for depth that is different  !
!    than the original coding.                                         !
!                                                                      !
!  * This version of the code has been modified by Bronwyn Cahill and  !
!    includes a semi-Lagrangian vertical sinking flux algorithm for    !
!    fecal material and bio_sediment subroutine which remineralizes    !
!    particulate nitrogen and returns it into dissolved nitrate pool.  !
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
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % SpecIr,                           &
     &                   FORCES(ng) % avcos,                            &
     &                   OCEAN(ng) % t)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!***********************************************************************
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w,                            &
     &                         SpecIr, avcos,                           &
     &                         t)
!***********************************************************************
!
      USE mod_param
      USE mod_biology
      USE mod_eclight
      USE mod_scalars
      USE mod_iounits
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
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: SpecIr(LBi:,LBj:,:)
      real(r8), intent(in) :: avcos(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: SpecIr(LBi:UBi,LBj:UBj,NBands)
      real(r8), intent(in) :: avcos(LBi:UBi,LBj:UBj,NBands)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Msink = 30

      integer :: Iter, Tindex, i, isink, ibio, id, itrc, j, k, ic, ks
      integer :: ibac, iband, idom, ifec, iphy, ipig
      integer :: Nsink

      integer, dimension(Msink) :: idsink
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), parameter :: MinVal = 0.0_r8

      real(r8) :: FV1, FV2, FV3, FV4, FV5, FV6, FV7, dtbio
      real(r8) :: DOC_lab, Ed_tot, Nup_max, aph442, aPHYN_wa
      real(r8) :: avgcos_min, par_b, par_bb, photo_DIC, photo_DOC
      real(r8) :: photo_decay, slope_AC, tChl, theta_m, total_photo
      real(r8) :: tot_ab, tot_b, tot_bb

      real(r8) :: Het_BAC
      real(r8) :: N_quota, RelDOC1, RelDON1, RelDOP1, RelFe
      real(r8) :: cff, cff1, cffL, cffR, cu, dltL, dltR

      real(r8), dimension(Msink) :: Wbio

      real(r8), dimension(4) :: Bac_G

      real(r8), dimension(NBands) :: dATT_sum

      real(r8), dimension(N(ng),NBands) :: specir_d
      real(r8), dimension(N(ng),NBands) :: avgcos, dATT

      real(r8), dimension(N(ng),Nphy) :: C2CHL, C2CHL_w
      real(r8), dimension(N(ng),Nphy) :: Gt_fl, Gt_ll, Gt_nl
      real(r8), dimension(N(ng),Nphy) :: Gt_sl, Gt_pl
      real(r8), dimension(N(ng),Nphy) :: alfa
      real(r8), dimension(N(ng),Nphy) :: pac_eff

      real(r8), dimension(N(ng),Nphy,Npig) :: Pigs_w

      integer, dimension(IminS:ImaxS) :: Keuphotic

      real(r8), dimension(IminS:ImaxS,N(ng)) :: E0_nz
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Ed_nz
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DOC_frac
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NitrBAC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NH4toNO3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NtoNBAC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NtoPBAC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NtoFeBAC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totDOC_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totDON_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totDOP_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totFe_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totNH4_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totNO3_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totPO4_d
      real(r8), dimension(IminS:ImaxS,N(ng)) :: totSiO_d

      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: GtBAC
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupDOC_ba
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupDON_ba
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupDOP_ba
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupFe_ba
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupNH4_ba
      real(r8), dimension(IminS:ImaxS,N(ng),Nbac) :: NupPO4_ba

      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: C2fALG
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: C2nALG
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: C2pALG
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: C2sALG
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: GtALG
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: GtALG_r
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupDOP
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupDON
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupFe
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupNH4
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupNO3
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupPO4
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: NupSiO
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: graz_act
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: mu_bar_f
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: mu_bar_n
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: mu_bar_p
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: mu_bar_s
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy) :: refuge

      real(r8), dimension(IminS:ImaxS,N(ng),Nfec) :: Regen_C
      real(r8), dimension(IminS:ImaxS,N(ng),Nfec) :: Regen_F
      real(r8), dimension(IminS:ImaxS,N(ng),Nfec) :: Regen_N
      real(r8), dimension(IminS:ImaxS,N(ng),Nfec) :: Regen_P
      real(r8), dimension(IminS:ImaxS,N(ng),Nfec) :: Regen_S

      real(r8), dimension(IminS:ImaxS,N(ng),NBands) :: specir_scal
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy,NBands) :: aPHYN_al
      real(r8), dimension(IminS:ImaxS,N(ng),Nphy,NBands) :: aPHYN_at

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_new

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
!
!=======================================================================
!  Add EcoSim Source/Sink terms.
!=======================================================================
!
!  Set internal time-stepping.
!
      dtbio=dt(ng)/REAL(BioIter(ng),r8)
!
!  Set vertical sinking identification and associated sinking velocity
!  arrays.
!
      ic=1
      DO ifec=1,Nfec
        idsink(ic)=iFecN(ifec)
        Wbio(ic)=WF(ifec,ng)
        ic=ic+1
        idsink(ic)=iFecC(ifec)
        Wbio(ic)=WF(ifec,ng)
        ic=ic+1
        idsink(ic)=iFecP(ifec)
        Wbio(ic)=WF(ifec,ng)
        ic=ic+1
        idsink(ic)=iFecS(ifec)
        Wbio(ic)=WF(ifec,ng)
        ic=ic+1
        idsink(ic)=iFecF(ifec)
        Wbio(ic)=WF(ifec,ng)
        ic=ic+1
      END DO

      DO iphy=1,Nphy
        idsink(ic)=iPhyN(iphy)
        Wbio(ic)=WS(iphy,ng)
        ic=ic+1
        idsink(ic)=iPhyC(iphy)
        Wbio(ic)=WS(iphy,ng)
        ic=ic+1
        idsink(ic)=iPhyP(iphy)
        Wbio(ic)=WS(iphy,ng)
        IF (iPhyS(iphy).ne.0) THEN
          ic=ic+1
          idsink(ic)=iPhyS(iphy)
          Wbio(ic)=WS(iphy,ng)
        END IF
        ic=ic+1
        idsink(ic)=iPhyF(iphy)
        Wbio(ic)=WS(iphy,ng)
        ic=ic+1
      END DO
      Nsink=ic-1
!
!-----------------------------------------------------------------------
!  Compute inverse thickness to avoid repeated divisions.
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!-----------------------------------------------------------------------
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO ibio=1,NBT
          itrc=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,itrc)=MAX(MinVal,t(i,j,k,nstp,itrc))
              Bio_old(i,k,itrc)=Bio(i,k,itrc)
!!
!! HGA - The new tendency terms were not initialized.  This gives
!!       unexpected behavior on different computers since a variable
!!       was used before it was assigned.  This may explain earlier
!!       problems with the algorithm.  Perhaps, this time-stepping can
!!       be modified latter to avoid unnecessary storage between
!!       Bio, Bio_old, and Bio_new.
!!
              Bio_new(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=t(i,j,k,nstp,isalt)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute temperature and salinity dependent variables.
!-----------------------------------------------------------------------
!
!  Refuge depth calculation.
!
        DO iphy=1,Nphy
          DO k=1,N(ng)
            DO i=Istr,Iend
              refuge(i,k,iphy)=MinRefuge(iphy,ng)
            END DO
          END DO
        END DO
!
!  Initialize fecal regeneration arrays (N, P, and Fe from Moore et al.,
!  DSRII 2001; silica is given by values from Bidle and Azam, Nature,
!  1999).
!
        IF (Regen_flag(ng)) THEN
          DO ifec=1,Nfec
            DO k=1,N(ng)
              DO i=Istr,Iend
                FV1=EXP(RegTfac(ifec,ng)*(Bio(i,k,itemp)-               &
     &                  RegTbase(ifec,ng)))
                Regen_C(i,k,ifec)=RegCR(ifec,ng)*FV1
                Regen_N(i,k,ifec)=RegNR(ifec,ng)*FV1
                Regen_P(i,k,ifec)=RegPR(ifec,ng)*FV1
                Regen_F(i,k,ifec)=RegFR(ifec,ng)*FV1
                Regen_S(i,k,ifec)=RegSR(ifec,ng)*FV1
              END DO
            END DO
          END DO
        END IF
!
!  Calculate temperature dependent growth rate.
!
        DO iphy=1,Nphy
          DO k=1,N(ng)
            DO i=Istr,Iend
              GtALG(i,k,iphy)=GtALG_max(iphy,ng)*                       &
     &                        EXP(PhyTfac(iphy,ng)*                     &
     &                            (Bio(i,k,itemp)-PhyTbase(iphy,ng)))
!
!  Calculate mu_bar for droop equation.
!
              FV1=maxC2nALG(iphy,ng)*(1.0_r8+GtALG(i,k,iphy))
              mu_bar_n(i,k,iphy)=GtALG(i,k,iphy)*                       &
     &                           FV1/(FV1-minC2nALG(iphy,ng))
              IF (HsSiO(iphy,ng).lt.LARGE) THEN
                FV1=maxC2SiALG(iphy,ng)*(1.0_r8+GtALG(i,k,iphy))
                mu_bar_s(i,k,iphy)=GtALG(i,k,iphy)*                     &
     &                             FV1/(FV1-minC2SiALG(iphy,ng))
              ELSE
                mu_bar_s(i,k,iphy)=LARGE
              END IF
              IF (HsPO4(iphy,ng).lt.LARGE) THEN
                FV1=maxC2pALG(iphy,ng)*(1.0_r8+GtALG(i,k,iphy))
                mu_bar_p(i,k,iphy)=GtALG(i,k,iphy)*                     &
     &                             FV1/(FV1-minC2pALG(iphy,ng))
              ELSE
                mu_bar_p(i,k,iphy)=LARGE
              END IF
              IF (HsFe(iphy,ng).lt.LARGE) THEN
                FV1=maxC2FeALG(iphy,ng)*(1.0_r8+GtALG(i,k,iphy))
                mu_bar_f(i,k,iphy)=GtALG(i,k,iphy)*                     &
     &                             FV1/(FV1-minC2FeALG(iphy,ng))
              ELSE
                mu_bar_f(i,k,iphy)=LARGE
              END IF
            END DO
          END DO
        END DO
!
!  Bacterial growth rate from Fasham et al., 1990.
!
        DO ibac=1,Nbac
          DO k=1,N(ng)
            DO i=Istr,Iend
              GtBAC(i,k,ibac)=GtBAC_max(ibac,ng)*                       &
     &                        EXP(BacTfac(ibac,ng)*                     &
     &                            (Bio(i,k,itemp)-BacTbase(ibac,ng)))
            END DO
          END DO
        END DO
!
!  Grazing rate calculation.
!  NOTE: ES1 included separation calculations for grazing beneath the
!        zone of refuge (250 m). This has been removed and may
!        result in differences in deeper waters.
!! Revisions, WPB 10/20/02.  New grazing formulation that is better
!! representation of basal loss rates and biomass accumulations.
!
        DO iphy=1,Nphy
          DO k=1,N(ng)
            DO i=Istr,Iend
              FV1=MAX(1.0_r8,(Bio(i,k,iPhyC(iphy))/refuge(i,k,iphy)))
              graz_act(i,k,iphy)=HsGRZ(iphy,ng)*LOG(FV1)
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Iterate biology source and sink terms.
!-----------------------------------------------------------------------
!
        ITER_LOOP : DO Iter=1,BioIter(ng)

          DO k=1,N(ng)
            DO i=Istr,Iend
              totNH4_d(i,k)=0.0_r8
              totNO3_d(i,k)=0.0_r8
              totPO4_d(i,k)=0.0_r8
              totSiO_d(i,k)=0.0_r8
              totFe_d (i,k)=0.0_r8
              totDOC_d(i,k)=0.0_r8
              totDON_d(i,k)=0.0_r8
              totDOP_d(i,k)=0.0_r8
            END DO
          END DO
          DO iphy=1,Nphy
            DO k=1,N(ng)
              DO i=Istr,Iend
                NupNH4(i,k,iphy)=0.0_r8
                NupNO3(i,k,iphy)=0.0_r8
                NupPO4(i,k,iphy)=0.0_r8
                NupSiO(i,k,iphy)=0.0_r8
                NupFe (i,k,iphy)=0.0_r8
                NupDON(i,k,iphy)=0.0_r8
                NupDOP(i,k,iphy)=0.0_r8
              END DO
            END DO
          END DO
!
!  Compute Ratio Arrays.
!  (Calculating only those that are accessed more than once.)
!
          DO iphy=1,Nphy
            DO k=1,N(ng)
              DO i=Istr,Iend
                C2nALG(i,k,iphy)=0.0_r8
                IF (Bio(i,k,iPhyN(iphy)).gt.0.0_r8) THEN
                  C2nALG(i,k,iphy)=Bio(i,k,iPhyC(iphy))/                &
     &                             Bio(i,k,iPhyN(iphy))
                END IF
                C2pALG(i,k,iphy)=0.0_r8
                IF (Bio(i,k,iPhyP(iphy)).gt.0.0_r8) THEN
                  C2pALG(i,k,iphy)=Bio(i,k,iPhyC(iphy))/                &
     &                             Bio(i,k,iPhyP(iphy))
                END IF
                C2sALG(i,k,iphy)=0.0_r8
                IF (iPhyS(iphy).gt.0) THEN
                  IF (Bio(i,k,iPhyS(iphy)).gt.0.0_r8) THEN
                    C2sALG(i,k,iphy)=Bio(i,k,iPhyC(iphy))/              &
     &                               Bio(i,k,iPhyS(iphy))
                  END IF
                END IF
                C2fALG(i,k,iphy)=0.0_r8
                IF (Bio(i,k,iPhyF(iphy)).gt.0.0_r8) THEN
                  C2fALG(i,k,iphy)=Bio(i,k,iPhyC(iphy))/                &
     &                             Bio(i,k,iPhyF(iphy))
                END IF
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Daylight Computations.
!-----------------------------------------------------------------------
!
! Initialize.
!
          DO i=Istr,Iend
            Ed_nz(i,N(ng))=0.0_r8
            E0_nz(i,N(ng))=0.0_r8
            Keuphotic(i)=N(ng)+1
            IF (SpecIr(i,j,21).gt.VSMALL) THEN
              DO k=1,N(ng)-1
                Ed_nz(i,k)=0.0_r8
                E0_nz(i,k)=0.0_r8
              END DO
              DO iband=1,NBands
                dATT_sum(iband)=0.0_r8
                DO k=1,N(ng)
                  dATT(k,iband)=0.0_r8
                END DO
                DO iphy=1,Nphy
                  DO k=1,N(ng)
                    aPHYN_at(i,k,iphy,iband)=0.0_r8
                    aPHYN_al(i,k,iphy,iband)=0.0_r8
                  END DO
                END DO
              END DO
!
!  Calculate average cosine zenith angle at surface.
!  (See equation 14 Morel, 1991 Prog. Ocean.)
!
              Ed_tot=0.0_r8
              DO iband=1,NBands
                Ed_tot=Ed_tot+SpecIr(i,j,iband)*DLAM
                avgcos(N(ng),iband)=avcos(i,j,iband)
              END DO
!
!  Total aph(442). adp(442) is set to 50% of aph(442).
!  NOTE: choosing sbands=9 which is band 442 using v8r16
!        sbands formulation. If spectral resolution changes, this
!        value must change!
!
              DO k=N(ng),1,-1
                IF (Ed_tot.ge.1.0_r8) THEN
                  aph442=0.0_r8
                  tChl=0.0_r8
                  DO iphy=1,Nphy
                    IF (Bio(i,k,iPhyC(iphy)).gt.0.0_r8) THEN
                      tChl=tChl+Bio(i,k,iPigs(iphy,ichl))
                      pac_eff(k,iphy)=1.0_r8
                      IF (b_PacEff(iphy,ng).gt.SMALL) THEN
                        FV2=Bio(i,k,iPigs(iphy,ichl))/                  &
     &                      (Bio(i,k,iPhyC(iphy))*12.0_r8)
                        pac_eff(k,iphy)=MAX(0.5_r8,                     &
     &                                      (MIN(1.0_r8,                &
     &                                           b_PacEff(iphy,ng)+     &
     &                                           mxPacEff(iphy,ng)*     &
     &                                           (FV2-                  &
     &                                            b_C2Cl(iphy,ng)))))
                      END IF
                      iband=9
                      DO ipig=1,Npig
                        IF (iPigs(iphy,ipig).gt.0) THEN
                          aph442=aph442+                                &
     &                           Bio(i,k,iPigs(iphy,ipig))*             &
     &                           apigs(ipig,iband)*pac_eff(k,iphy)
                        END IF
                      END DO
                    END IF
                  END DO
!
!  Calculate absorption.
!  Calculating phytoplankton absorption for attentuation calculation.
!  NOTE: 12 factor to convert to ugrams (mg m-3)
!
                  aph442=0.5_r8*aph442
                  DO iband=1,NBands
                    tot_ab=0.0_r8
                    DO iphy=1,Nphy
                      DO ipig=1,Npig
                        IF (iPigs(iphy,ipig).gt.0) THEN
                          aPHYN_at(i,k,iphy,iband)=                     &
     &                                      aPHYN_at(i,k,iphy,iband)+   &
     &                                      Bio(i,k,iPigs(iphy,ipig))*  &
     &                                      apigs(ipig,iband)*          &
     &                                      pac_eff(k,iphy)
                        END IF
                      END DO
                      tot_ab=tot_ab+aPHYN_at(i,k,iphy,iband)
!
!  Removing absorption due to PPC for "alfa" calculation.
!
                      ipig=5
                      IF (iPigs(iphy,ipig).gt.0) THEN
                        aPHYN_al(i,k,iphy,iband)=                       &
     &                                    aPHYN_at(i,k,iphy,iband)-     &
     &                                    Bio(i,k,iPigs(iphy,ipig))*    &
     &                                    apigs(ipig,iband)*            &
     &                                    pac_eff(k,iphy)
                      END IF
                    END DO
!
!  Adding detrital absorption.
!
                    tot_ab=tot_ab+                                      &
     &                     aph442*EXP(0.011_r8*(442.0_r8-               &
     &                                (397.0_r8+REAL(iband,r8)*DLAM)))
!
!  Calculate CDOC absorption.
!  NOTE: 12 factor is to convert ugrams per liter, and 0.001 converts
!        to mg/liter.  Specific absorption
!        coefficients were calculated as m-1 / (mg DOC/liters sw).
!        net factor = (12*0.001) = 0.012
!
                    tot_ab=tot_ab+                                      &
     &                     0.012_r8*(Bio(i,k,iCDMC(ilab))*              &
     &                               aDOC(ilab,iband)+                  &
     &                               Bio(i,k,iCDMC(irct))*              &
     &                               aDOC(irct,iband))+                 &
     &                     awater(iband)
!
!  Calculate scattering and backscattering (see equation 19 Morel, 1991,
!  Prog. Ocean). Morel, 1988 puts spectral dependency in backscattering.
!  Since Morel (1991) does not have a backscattering equation, use 1988
!  paper. Morel 2001 has slight adjustment 0.01, rather than 0.02.
!  This was altered, but never tested in ROMS 1.8 on 03/08/03.
!
                    par_b =0.3_r8*(tChl**0.62_r8)
                    par_bb=0.0_r8
                    IF (tChl.gt.0.0_r8) THEN
                      par_bb=par_b*(0.002_r8+0.01_r8*                   &
     &                              (0.5_r8-0.25_r8*LOG10(tChl))*       &
     &                              wavedp(iband))
                    END IF
                    par_bb=MAX(par_bb,0.0_r8)
!
!  However, for omega0 calculation, par_b must be spectral, so use
!  dependency from Sathy and Platt 1988
!
                    tot_b=bwater(iband)+par_b*wavedp(iband)
!
!  Morel, 1988 instead of 1991. See methods
!
                    tot_bb=0.5_r8*bwater(iband)+par_bb
!
!  Sathy and Platt JGR 1988.  This is set with the average cosine of
!  the box above, and used to calculate a new avgcos for this level.
!  This new average cosine is then used to recalculate the attenuation
!  coefficient
!
                    dATT(k,iband)=(tot_ab+tot_bb)/avgcos(k,iband)
!
!  See Mobley, 1995 for graphical depiction of this equation.
!
                    avgcos_min=avgcos(k,iband)+                         &
     &                         (0.5_r8-avgcos(k,iband))*                &
     &                         (tot_b/(tot_ab+tot_b))
!
!  Calculate average cosine. Linear fit to average cosine versus optical
!  depth relationship. The FV1 calculation keeps the denominator of the
!  slope calculation from going negative and above 1.
!
                    FV1=MAX(1.0_r8,                                     &
     &                      7.0_r8-dATT(k,iband)*ABS(z_r(i,j,k)))
                    slope_AC =MIN(0.0_r8,                               &
     &                            (avgcos_min-avgcos(k,iband))/FV1)
                    avgcos(k,iband)=avgcos(k,iband)+                    &
     &                             slope_AC*dATT(k,iband)*Hz(i,j,k)
                    dATT(k,iband)=(tot_ab+tot_bb)/avgcos(k,iband)
!
!  Set avgcos for next level.
!
                    IF (k.ne.1) THEN
                      avgcos(k-1,iband)=avgcos(k,iband)
                    END IF
!
!  Calculate spectral irradiance with depth.
!
                    FV1=dATT(k,iband)*Hz(i,j,k)
                    FV2=dATT_sum(iband)+0.5_r8*FV1
                    dATT_sum(iband)=dATT_sum(iband)+FV1
                    specir_d(k,iband)=SpecIr(i,j,iband)*                &
     &                                EXP(-FV2)*DLAM
!
!  Calculate spectral scalar irradiance.  Morel, 1991 Prog. Ocean.
!
                    specir_scal(i,k,iband)=specir_d(k,iband)*           &
     &                                     (dATT(k,iband)/tot_ab)
                    E0_nz(i,k)=E0_nz(i,k)+specir_scal(i,k,iband)
!
!  Calculate Ed_nz.
!
                    Ed_nz(i,k)=Ed_nz(i,k)+specir_d(k,iband)
                  END DO
                  Ed_tot=E0_nz(i,k)
!
!  Set bottom of the euphotic zone.
!
                  Keuphotic(i)=k
                END IF
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Bacterial nutrient uptake.
!-----------------------------------------------------------------------
!
          DO ibac=1,Nbac
            DO k=1,N(ng)
              DO i=Istr,Iend
!
!  DOM uptake.
!
                IF ((Bio(i,k,iDOMC(ilab)).gt.0.0_r8).and.               &
     &              (Bio(i,k,iDOMN(ilab)).gt.0.0_r8).and.               &
     &              (Bio(i,k,iDOMP(ilab)).gt.0.0_r8)) THEN
                  NupDOC_ba(i,k,ibac)=GtBAC(i,k,ibac)*                  &
     &                                Bio(i,k,iBacC(ibac))*             &
     &                                I_Bac_Ceff(ng)*                   &
     &                                (Bio(i,k,iDOMC(ilab))/            &
     &                                (HsDOC_ba(ibac,ng)+               &
     &                                 Bio(i,k,iDOMC(ilab))))
                  NupDON_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)*              &
     &                                Bio(i,k,iDOMN(ilab))/             &
     &                                Bio(i,k,iDOMC(ilab))
                  NupDOP_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)*              &
     &                                Bio(i,k,iDOMP(ilab))/             &
     &                                Bio(i,k,iDOMC(ilab))
                ELSE
                  NupDOC_ba(i,k,ibac)=0.0_r8
                  NupDON_ba(i,k,ibac)=0.0_r8
                  NupDOP_ba(i,k,ibac)=0.0_r8
                END IF
                totDOC_d(i,k)=totDOC_d(i,k)+NupDOC_ba(i,k,ibac)
                totDON_d(i,k)=totDON_d(i,k)+NupDON_ba(i,k,ibac)
                totDOP_d(i,k)=totDOP_d(i,k)+NupDOP_ba(i,k,ibac)
!
!  NH4 uptake.
!
                NupNH4_ba(i,k,ibac)=GtBAC(i,k,ibac)*                    &
     &                              Bio(i,k,iBacN(ibac))*               &
     &                              Bio(i,k,iNH4_)/                     &
     &                              (HsNH4_ba(ibac,ng)+Bio(i,k,iNH4_))
                totNH4_d(i,k)=totNH4_d(i,k)+NupNH4_ba(i,k,ibac)
!
!  PO4 uptake.
!
                NupPO4_ba(i,k,ibac)=GtBAC(i,k,ibac)*                    &
     &                              Bio(i,k,iBacP(ibac))*               &
     &                              Bio(i,k,iPO4_)/                     &
     &                              (HsPO4_ba(ibac,ng)+Bio(i,k,iPO4_))
                totPO4_d(i,k)=totPO4_d(i,k)+NupPO4_ba(i,k,ibac)
!
!  Fe uptake.
!
                NupFe_ba(i,k,ibac)=GtBAC(i,k,ibac)*                     &
     &                             Bio(i,k,iBacF(ibac))*                &
     &                             Bio(i,k,iFeO_)/                      &
     &                             (HsFe_ba(ibac,ng)+Bio(i,k,iFeO_))
                totFe_d(i,k)=totFe_d(i,k)+NupFe_ba(i,k,ibac)
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton dark nutrient uptake.
!-----------------------------------------------------------------------
!
          DO iphy=1,Nphy
            DO k=1,N(ng)
              DO i=Istr,Iend
                IF (C2nALG(i,k,iphy).gt.C2nALGminABS(iphy,ng)) THEN
!
!  NOTE: these are being saved to test for total nutrient uptake.
!        If nutrient uptake is greater than maximum nutrient, then
!        each of the uptakes are reduced by their fractional contri-
!        bution to the total.
!
                  Nup_max=GtALG(i,k,iphy)
                  NupNO3(i,k,iphy)=(Bio(i,k,iNO3_)/                     &
     &                             (HsNO3(iphy,ng)+Bio(i,k,iNO3_))*     &
     &                             EXP(-BET_(iphy,ng)*Bio(i,k,iNH4_)))
                  NupNH4(i,k,iphy)=Bio(i,k,iNH4_)/                      &
     &                             (HsNH4(iphy,ng)+Bio(i,k,iNH4_))
!
!  Test that Wroblewski equation does not exceed 1.0.
!
                  FV1=NupNO3(i,k,iphy)+NupNH4(i,k,iphy)
                  IF (FV1.gt.1.0_r8) THEN
                    FV1=1.0_r8/FV1
                    NupNO3(i,k,iphy)=NupNO3(i,k,iphy)*FV1
                    NupNH4(i,k,iphy)=NupNH4(i,k,iphy)*FV1
                  END IF
!
!  Change from percentage of maximum to mass per second.
!
                  FV1=Nup_max*Bio(i,k,iPhyN(iphy))
                  NupNO3(i,k,iphy)=NupNO3(i,k,iphy)*FV1
                  NupNH4(i,k,iphy)=NupNH4(i,k,iphy)*FV1
!
!  Test for DON uptake.
!
                  IF (C2nALG(i,k,iphy).gt.C2nNupDON(iphy,ng)) THEN
                    NupDON(i,k,iphy)=FV1*                               &
     &                               Bio(i,k,iDOMN(ilab))/              &
     &                               (HsDON(iphy,ng)+                   &
     &                                Bio(i,k,iDOMN(ilab)))
                  END IF
!
!  Accumulate total demand for nutrients.
!
                  totNO3_d(i,k)=totNO3_d(i,k)+NupNO3(i,k,iphy)
                  totNH4_d(i,k)=totNH4_d(i,k)+NupNH4(i,k,iphy)
                  totDON_d(i,k)=totDON_d(i,k)+NupDON(i,k,iphy)
                END IF
!
!  Dark silica uptake, min C2Si test.
!  The LARGE test can be removed after testing phase.
!
                IF (HsSiO(iphy,ng).lt.LARGE) THEN
                  IF (C2sALG(i,k,iphy).gt.C2SiALGminABS(iphy,ng)) THEN
                    Nup_max=GtALG(i,k,iphy)
                    NupSiO(i,k,iphy)=Bio(i,k,iSiO_)/                    &
     &                               (HsSiO(iphy,ng)+Bio(i,k,iSiO_))
!
!  Change from percentage of maximum to mass per second.
!
                    IF (iPhyS(iphy).gt.0) THEN
                      FV1=Nup_max*Bio(i,k,iPhyS(iphy))
                      NupSiO(i,k,iphy)=NupSiO(i,k,iphy)*FV1
                    ELSE
                      NupSiO(i,k,iphy)=0.0_r8
                    END IF
!
!  Accumulate total demand for nutrients.
!
                    totSiO_d(i,k)=totSiO_d(i,k)+NupSiO(i,k,iphy)
                  END IF
                END IF
!
!  Dark phophorus uptake, min C2P test.
!  The LARGE test can be removed after testing phase.
!
                IF (HsPO4(iphy,ng).lt.LARGE) THEN
                  IF (C2pALG(i,k,iphy).gt.C2pALGminABS(iphy,ng)) THEN
                    Nup_max=GtALG(i,k,iphy)
                    NupPO4(i,k,iphy)=Bio(i,k,iPO4_)/                    &
     &                               (HsPO4(iphy,ng)+Bio(i,k,iPO4_))
!
!  Change from percentage of maximum to mass per second.
!
                    FV1=Nup_max*Bio(i,k,iPhyP(iphy))
                    NupPO4(i,k,iphy)=NupPO4(i,k,iphy)*FV1
!
!  Test for alk. phosphatase
!
                    IF (C2pALG(i,k,iphy).gt.C2pALKPHOS(iphy,ng)) THEN
                      NupDOP(i,k,iphy)=FV1*                             &
                                       Bio(i,k,iDOMP(ilab))/            &
     &                                 (HsDOP(iphy,ng)+                 &
     &                                  Bio(i,k,iDOMP(ilab)))
                    END IF
!
!  Accumulate total demand for nutrients.
!
                    totPO4_d(i,k)=totPO4_d(i,k)+NupPO4(i,k,iphy)
                    totDOP_d(i,k)=totDOP_d(i,k)+NupDOP(i,k,iphy)
                  END IF
                END IF
!
!  Dark iron uptake, min C2Fe test.
!  The LARGE test can be removed after testing phase.
!
                IF (HsFe(iphy,ng).lt.LARGE) THEN
                  IF (C2fALG(i,k,iphy).gt.C2FeALGminABS(iphy,ng)) THEN
                    Nup_max=GtALG(i,k,iphy)
                    NupFe(i,k,iphy)=Bio(i,k,iFeO_)/                     &
     &                              (HsFe(iphy,ng)+Bio(i,k,iFeO_))
!
!  Change from percentage of maximum to mass per second.
!
                    FV1=Nup_max*Bio(i,k,iPhyF(iphy))
                    NupFe(i,k,iphy)=NupFe(i,k,iphy)*FV1
!
!  Accumulate total demand for nutrients.
!
                    totFe_d(i,k)=totFe_d(i,k)+NupFe(i,k,iphy)
                  END IF
                END IF
              END DO
            END DO
          END DO
!
!  Calculate bacterial nitrification as a Michaelis-Menton function
!  of ambient NH4 concentration, beneath the euphotic zone (light
!  inhibits nitrification).
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              NitrBAC(i,k)=0.0_r8
              NH4toNO3(i,k)=0.0_r8
              NtoNBAC(i,k)=0.0_r8
              NtoPBAC(i,k)=0.0_r8
              NtoFeBAC(i,k)=0.0_r8
              IF (k.lt.Keuphotic(i)) THEN
                NH4toNO3(i,k)=RtNIT(ng)*                                &
     &                        Bio(i,k,iNH4_)/(HsNIT(ng)+Bio(i,k,iNH4_))
!
!  Nitrification fixes DIC into POC.
!  Conversion factor of 7.0 from Kaplan 1983 "Nitrogen in the Sea"
!  factor equals (1.0 / (7.0 * C2nBAC)). Adds NH4 uptake as biomass.
!
                NitrBAC(i,k)=NH4toNO3(i,k)/7.0_r8
                NtoNBAC(i,k)=NitrBAC(i,k)*N2cBAC(ng)
                NtoPBAC(i,k)=NitrBAC(i,k)*P2cBAC(ng)
                NtoFeBAC(i,k)=NitrBAC(i,k)*Fe2cBAC(ng)
                totNH4_d(i,k)=totNH4_d(i,k)+NH4toNO3(i,k)+NtoNBAC(i,k)
                totPO4_d(i,k)=totPO4_d(i,k)+NtoPBAC(i,k)
                totFe_d (i,k)=totFe_d (i,k)+NtoFeBAC(i,k)
              END IF
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Test that total nutrient demand does not exceed supply.  If it does
!  total demand is normalized to the total supply. Each species demand
!  is reduced to its weighted average percentage of the supply.
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              FV2=totNO3_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iNO3_)) THEN
                FV1=(Bio(i,k,iNO3_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupNO3(i,k,iphy)=NupNO3(i,k,iphy)*FV1
                END DO
              END IF
!
              FV2=totNH4_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iNH4_)) THEN
                FV1=(Bio(i,k,iNH4_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupNH4(i,k,iphy)=NupNH4(i,k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupNH4_ba(i,k,ibac)=NupNH4_ba(i,k,ibac)*FV1
                END DO
                NH4toNO3(i,k)=NH4toNO3(i,k)*FV1
                NtoNBAC(i,k)=NtoNBAC(i,k)*FV1
              END IF
!
              FV2=totSiO_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iSiO_)) THEN
                FV1=(Bio(i,k,iSiO_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupSiO(i,k,iphy)=NupSiO(i,k,iphy)*FV1
                END DO
              END IF
!
              FV2=totPO4_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iPO4_)) THEN
                FV1=(Bio(i,k,iPO4_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupPO4(i,k,iphy)=NupPO4(i,k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupPO4_ba(i,k,ibac)=NupPO4_ba(i,k,ibac)*FV1
                END DO
                NtoPBAC(i,k)=NtoPBAC(i,k)*FV1
              END IF
!
              FV2=totFe_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iFeO_)) THEN
                FV1=(Bio(i,k,iFeO_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupFe(i,k,iphy)=NupFe(i,k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupFe_ba(i,k,ibac)=NupFe_ba(i,k,ibac)*FV1
                END DO
                NtoFeBAC(i,k)=NtoFeBAC(i,k)*FV1
              END IF
!
!  Bacteria are the only group to take up DOC.  Remove BAC DON and
!  BAC DOP uptake from total uptake; adjust uptake and add back.
!
              FV2=totDOC_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iDOMC(ilab))) THEN
                FV1=(Bio(i,k,iDOMC(ilab))-VSMALL)/FV2
                totDOC_d(i,k)=totDOC_d(i,k)*FV1
                DO ibac=1,Nbac
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)*FV1
                  totDON_d(i,k)=totDON_d(i,k)-NupDON_ba(i,k,ibac)
                  NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)*FV1
                  totDON_d(i,k)=totDON_d(i,k)+NupDON_ba(i,k,ibac)
                  totDOP_d(i,k)=totDOP_d(i,k)-NupDOP_ba(i,k,ibac)
                  NupDOP_ba(i,k,ibac)=NupDOP_ba(i,k,ibac)*FV1
                  totDOP_d(i,k)=totDOP_d(i,k)+NupDOP_ba(i,k,ibac)
                END DO
              END IF
!
!  Remove BAC DON uptake from total uptake; adjust uptake and add back.
!
              FV2=totDON_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iDOMN(ilab))) THEN
                FV1=(Bio(i,k,iDOMN(ilab))-VSMALL)/FV2
                totDON_d(i,k)=totDON_d(i,k)*FV1
                totDOC_d(i,k)=totDOC_d(i,k)*FV1
                DO iphy=1,Nphy
                  NupDON(i,k,iphy)=NupDON(i,k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)*FV1
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)*FV1
                  totDOP_d(i,k)=totDOP_d(i,k)-NupDOP_ba(i,k,ibac)
                  NupDOP_ba(i,k,ibac)=NupDOP_ba(i,k,ibac)*FV1
                  totDOP_d(i,k)=totDOP_d(i,k)+NupDOP_ba(i,k,ibac)
                END DO
              END IF
!
!  Remove BAC DOP uptake from total uptake; adjust uptake and add back.
!
              FV2=totDOP_d(i,k)*dtbio
              IF (FV2.gt.Bio(i,k,iDOMP(ilab))) THEN
                FV1=(Bio(i,k,iDOMP(ilab))-VSMALL)/FV2
                totDOP_d(i,k)=totDOP_d(i,k)*FV1
                totDOC_d(i,k)=totDOC_d(i,k)*FV1
                DO iphy=1,Nphy
                  NupDOP(i,k,iphy)=NupDOP(i,k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupDOP_ba(i,k,ibac)=NupDOP_ba(i,k,ibac)*FV1
                  totDON_d(i,k)=totDON_d(i,k)-NupDON_ba(i,k,ibac)
                  NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)*FV1
                  totDON_d(i,k)=totDON_d(i,k)+NupDON_ba(i,k,ibac)
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)*FV1
                END DO
              END IF
            END DO
          END DO
!
!  Increase particulate nutrients by the amount of the uptake.
!
          DO iphy=1,Nphy
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio_new(i,k,iPhyN(iphy))=Bio_new(i,k,iPhyN(iphy))+      &
     &                                   NupNO3(i,k,iphy)+              &
     &                                   NupNH4(i,k,iphy)+              &
     &                                   NupDON(i,k,iphy)
                Bio_new(i,k,iPhyP(iphy))=Bio_new(i,k,iPhyP(iphy))+      &
     &                                   NupPO4(i,k,iphy)+              &
     &                                   NupDOP(i,k,iphy)
                Bio_new(i,k,iPhyF(iphy))=Bio_new(i,k,iPhyF(iphy))+      &
     &                                   NupFe(i,k,iphy)
                IF (iPhyS(iphy).gt.0) THEN
                  Bio_new(i,k,iPhyS(iphy))=Bio_new(i,k,iPhyS(iphy))+    &
     &                                     NupSiO(i,k,iphy)
                END IF
!
!  Update nutrient arrays for growth and budgets. Bacterial uptake
!  included below.
!
                Bio_new(i,k,iNO3_)=Bio_new(i,k,iNO3_)-                  &
     &                             NupNO3(i,k,iphy)
                Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)-                  &
     &                             NupNH4(i,k,iphy)
                Bio_new(i,k,iSiO_)=Bio_new(i,k,iSiO_)-                  &
     &                             NupSiO(i,k,iphy)
                Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)-                  &
     &                             NupPO4(i,k,iphy)
                Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)-                  &
     &                             NupFe (i,k,iphy)
                Bio_new(i,k,iDOMN(ilab))=Bio_new(i,k,iDOMN(ilab))-      &
     &                                   NupDON(i,k,iphy)
                Bio_new(i,k,iDOMP(ilab))=Bio_new(i,k,iDOMP(ilab))-      &
     &                                   NupDOP(i,k,iphy)
              END DO
            END DO
          END DO
!
!  Nitrification fixes DIC into DOC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)-                    &
     &                           NitrBAC(i,k)
            END DO
          END DO
!
!  Add nitrifying bacteria biomass to heterotrophic bacteria biomass.
!  Adding PON, POP, POFe to BacC arrays at current C2_BAC ratios.
!
          DO ibac=1,Nbac
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio_new(i,k,iBacC(ibac))=Bio_new(i,k,iBacC(ibac))+      &
     &                                   NitrBAC(i,k)
                Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))+      &
     &                                   NtoNBAC(i,k)
                Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))+      &
     &                                   NtoPBAC(i,k)
                Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))+      &
     &                                   NtoFeBAC(i,k)
              END DO
            END DO
          END DO
!
!  Update nutrient arrays for nitrification.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_new(i,k,iNO3_)=Bio_new(i,k,iNO3_)+                    &
     &                           NH4toNO3(i,k)
              Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)-                    &
     &                           (NH4toNO3(i,k)+NtoNBAC(i,k))
              Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)-                    &
     &                           NtoPBAC(i,k)
              Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)-                    &
     &                           NtoFeBAC(i,k)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Light mediated carbon growth.
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
            DO k=N(ng),Keuphotic(i),-1
              DO iphy=1,Nphy
                IF (Bio(i,k,iPhyC(iphy)).gt.0.0_r8) THEN
!
!  Calculate weighted average spectral absorption.
!
                  aPHYN_wa=0.0_r8
                  DO iband=1,NBands
                    aPHYN_wa=aPHYN_wa+(aPHYN_al(i,k,iphy,iband)*        &
     &                                 specir_scal(i,k,iband))
                  END DO
!
!  If Keuphotic(i) < N+1, and E0_nz(i,k)=0, this will cause pigments to
!  blow up. This should never happen, unless Keuphotic is not calcuated
!  properly. WPB
!
                  aPHYN_wa=aPHYN_wa/E0_nz(i,k)
!
!  Calculate "alfa" for HTAN function of P vs. I.
!  (conversion:  Ein/microEin * 10e3)
!
                  alfa(k,iphy)=(aPHYN_wa/Bio(i,k,iPhyC(iphy)))*         &
     &                          qu_yld(iphy,ng)*0.001_r8
!
!  Light limited growth rate.
!
                  FV1=MAX(0.0_r8,E0_nz(i,k)-E0_comp(iphy,ng))
                  FV2=E0_nz(i,k)-E0_inhib(iphy,ng)
                  IF (FV2.gt.0.0_r8) THEN
                    Gt_ll(k,iphy)=GtALG(i,k,iphy)*                      &
     &                            TANH(alfa(k,iphy)*FV1/                &
     &                            GtALG(i,k,iphy))*                     &
     &                            EXP(-inhib_fac(iphy,ng)*FV2)
                  ELSE
                    Gt_ll(k,iphy)=GtALG(i,k,iphy)*                      &
     &                            TANH(alfa(k,iphy)*FV1/                &
     &                            GtALG(i,k,iphy))
                  END IF
!
!  Nutrient limited growth rates.
!
!  REMEMBER that sinking speed to be set by gradient of limiting
!       nutrient, allowing for negative sinking. Try storing growth
!       rate terms in an array and using MAXLOC for if test.
!
!  Nitrogen limited growth rate.
!
                  IF (Bio(i,k,iPhyN(iphy)).gt.0.0_r8) THEN
                    FV1=Bio(i,k,iPhyC(iphy))/                           &
     &                  (Bio(i,k,iPhyN(iphy))+Bio_new(i,k,iPhyN(iphy)))
                    Gt_nl(k,iphy)=mu_bar_n(i,k,iphy)*                   &
     &                            (1.0_r8-ImaxC2nALG(iphy,ng)*FV1)
                    Gt_nl(k,iphy)=MAX(0.0_r8,                           &
     &                                MIN(Gt_nl(k,iphy),                &
     &                                    GtALG(i,k,iphy)))
                  END IF
!
!  Silica limited growth rate.
!  Testing for silica incorporation.
!
                  IF (iPhyS(iphy).gt.0) THEN
                    IF ((HsSiO(iphy,ng).lt.LARGE).and.                  &
     &                  (Bio(i,k,iPhyS(iphy)).gt.0.0_r8)) THEN
                      FV1=Bio(i,k,iPhyC(iphy))/                         &
     &                    (Bio(i,k,iPhyS(iphy))+                        &
     &                     Bio_new(i,k,iPhyS(iphy)))
                      Gt_sl(k,iphy)=mu_bar_s(i,k,iphy)*                 &
     &                              (1.0_r8-ImaxC2SiALG(iphy,ng)*FV1)
                      Gt_sl(k,iphy)=MAX(0.0_r8,                         &
     &                                  MIN(Gt_sl(k,iphy),              &
     &                                      GtALG(i,k,iphy)))
                    ELSE
                      Gt_sl(k,iphy)=LARGE
                    END IF
                  ELSE
                    Gt_sl(k,iphy)=LARGE
                  END IF
!
!  Phosphorus limited growth rate.
!
                  IF ((HsPO4(iphy,ng).lt.LARGE).and.                    &
     &                (Bio(i,k,iPhyP(iphy)).gt.0.0_r8)) THEN
                    FV1=Bio(i,k,iPhyC(iphy))/                           &
     &                  (Bio(i,k,iPhyP(iphy))+Bio_new(i,k,iPhyP(iphy)))
                    Gt_pl(k,iphy)=mu_bar_p(i,k,iphy)*                   &
     &                            (1.0_r8-ImaxC2pALG(iphy,ng)*FV1)
                    Gt_pl(k,iphy)=MAX(0.0_r8,                           &
     &                                MIN(Gt_pl(k,iphy),                &
     &                                    GtALG(i,k,iphy)))
                  ELSE
                    Gt_pl(k,iphy)=LARGE
                  END IF
!
!  Iron limited growth rate
!
                  IF ((HsFe(iphy,ng).lt.LARGE).and.                     &
     &                (Bio(i,k,iPhyF(iphy)).gt.0.0_r8)) THEN
                    FV1=Bio(i,k,iPhyC(iphy))/                           &
     &                  (Bio(i,k,iPhyF(iphy))+Bio_new(i,k,iPhyF(iphy)))
                    Gt_fl(k,iphy)=mu_bar_f(i,k,iphy)*                   &
     &                            (1.0_r8-ImaxC2FeALG(iphy,ng)*FV1)
                    Gt_fl(k,iphy)=MAX(0.0_r8,                           &
     &                                MIN(Gt_fl(k,iphy),                &
     &                                    GtALG(i,k,iphy)))
                  ELSE
                    Gt_fl(k,iphy)=LARGE
                  END IF
!
!  Realized growth rate is minimum of light or nutrient limited rate.
!
                  GtALG_r(i,k,iphy)=MIN(Gt_ll(k,iphy),Gt_nl(k,iphy),    &
     &                                  Gt_sl(k,iphy),Gt_pl(k,iphy),    &
     &                                  Gt_fl(k,iphy))
                  IF (GtALG_r(i,k,iphy).ge.LARGE) THEN
                    GtALG_r(i,k,iphy)=0.0_r8
                  END IF
!
!  Carbon growth calculations.
!
                  FV1=Bio(i,k,iPhyC(iphy))*GtALG_r(i,k,iphy)
                  Bio_new(i,k,iPhyC(iphy))=Bio_new(i,k,iPhyC(iphy))+    &
     &                                     FV1
                  Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)-                &
     &                               FV1
!
!  Pigment growth calculations.
!
                  DO ipig=1,Npig
                    IF (iPigs(iphy,ipig).gt.0) THEN
                      itrc=iPigs(iphy,ipig)
                      IF (Bio(i,k,iPhyC(iphy)).gt.0.0_r8) THEN
                        FV1=Bio(i,k,itrc)*GtALG_r(i,k,iphy)
                        Bio_new(i,k,itrc)=Bio_new(i,k,itrc)+FV1
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Bacterioplankton carbon growth terms.
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Het_BAC=0.0_r8
              RelDOC1=0.0_r8
              RelDON1=0.0_r8
              RelDOP1=0.0_r8
              RelFe=0.0_r8
!
!  NOTE: Only DOC2/DON2 formation is in this section.
!        Take colored excretion off the top. 03/18/00
!        also, not excreting any DOP or Fe
!  REMEMBER, if excreting DOP and Fe, must address changes in growth if
!        tests. (see DON equations). 03/21/00.
!
              DO ibac=1,Nbac
                FV1=NupDOC_ba(i,k,ibac)*ExBAC_c(ng)*                    &
     &              (1.0_r8-cDOCfrac_c(irct,ng))
                FV2=NupDOC_ba(i,k,ibac)*ExBAC_c(ng)*                    &
     &              cDOCfrac_c(irct,ng)
                FV3=NupDON_ba(i,k,ibac)*ExBAC_n(ng)
!
                Bio_new(i,k,iDOMC(irct))=Bio_new(i,k,iDOMC(irct))+      &
     &                                   FV1
                Bio_new(i,k,iCDMC(irct))=Bio_new(i,k,iCDMC(irct))+      &
     &                                   FV2
                Bio_new(i,k,iDOMN(irct))=Bio_new(i,k,iDOMN(irct))+      &
     &                                   FV3
!
!  As we are taking it off the top, must remove from DOMN1 now. No other
!  organisms use DOMC1, so net term (totDOC_d) can be used in budgeting
!  below. This saves cycles, but makes code difficult to read. WPB
!
                Bio_new(i,k,iDOMN(ilab))=Bio_new(i,k,iDOMN(ilab))-      &
     &                                   FV3
!
!  Remove from uptake.
!
                NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)-                &
     &                              (FV1+FV2)
                NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)-                &
     &                              FV3
!
!  Determine growth limitation. Assuming 100% efficiency for N, P, Fe.
!  If DOMC=0, or DOMN=0, or DOMP=0, then NupDOC_ba = NupDON_ba =
!  NupDOP_ba = 0 and none of the divisions below are accessed. WPB
!
                Bac_G(1)=NupDOC_ba(i,k,ibac)*Bac_Ceff(ng)
                Bac_G(2)=(NupDON_ba(i,k,ibac)+                          &
     &                    NupNH4_ba(i,k,ibac))*                         &
     &                   C2nBAC(ng)
                Bac_G(3)=(NupDOP_ba(i,k,ibac)+                          &
     &                    NupPO4_ba(i,k,ibac))*                         &
     &                   C2pBAC(ng)
                Bac_G(4)=NupFe_ba(i,k,ibac)*C2FeBAC(ng)
!
!  Energy limited case. All excess nutrients returned in inorganic form.
!
                IF ((Bac_G(1).le.Bac_G(2)).and.                         &
     &              (Bac_G(1).le.Bac_G(3)).and.                         &
     &              (Bac_G(1).le.Bac_G(4))) THEN
                  Het_BAC=Bac_G(1)
                  FV1=Bac_G(1)*N2cBAC(ng)
                  FV2=Bac_G(1)*P2cBAC(ng)
                  FV3=Bac_G(1)*Fe2cBAC(ng)
                  Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))+    &
     &                                     FV1
                  Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))+    &
     &                                     FV2
                  Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays should probably now be negative. If NH4 or PO4 is
!  positive, then there is some uptake of inorganic forms, but this
!  value will be less than the original Nup value because of IF test.
!
                  NupNH4_ba(i,k,ibac)=FV1-NupDON_ba(i,k,ibac)
                  NupPO4_ba(i,k,ibac)=FV2-NupDOP_ba(i,k,ibac)
!
!  Because Fe is considered to be all inorganic, only net uptake of Fe
!  is needed.
!
                  RelFe=NupFe_ba(i,k,ibac)-FV3
                  NupFe_ba(i,k,ibac)=FV3
!
!  Nitrogen limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE IF ((Bac_G(2).le.Bac_G(3)).and.                    &
     &                   (Bac_G(2).le.Bac_G(4))) THEN
                  Het_BAC=Bac_G(2)
                  FV2=Bac_G(2)*P2cBAC(ng)
                  FV3=Bac_G(2)*Fe2cBAC(ng)
                  Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))+    &
     &                                     (NupDON_ba(i,k,ibac)+        &
     &                                      NupNH4_ba(i,k,ibac))
                  Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))+    &
     &                                     FV2
                  Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(2))*I_Bac_Ceff(ng)
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DOP from C2pDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c(ng)*                                &
!!   &                DOC_frac(i,k)*                                    &
     &                Bio(i,k,iDOMP(ilab))/                             &
     &                Bio(i,k,iDOMC(ilab))
                  FV5=FV2-(NupDOP_ba(i,k,ibac)+                         &
                           NupPO4_ba(i,k,ibac)-FV4)
!
!  If FV5 is positive then released DOP is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDOP1=FV4
                    NupPO4_ba(i,k,ibac)=NupPO4_ba(i,k,ibac)+FV5
                  ELSE
                    RelDOP1=FV4-FV5
                  END IF
                  NupDOP_ba(i,k,ibac)=NupDOP_ba(i,k,ibac)-RelDOP1
!
!  Release Fe.
!
                  RelFe=NupFe_ba(i,k,ibac)-FV3
                  NupFe_ba(i,k,ibac)=FV3
!
!  Phosphorous limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE IF (Bac_G(3).le.Bac_G(4)) THEN
                  Het_BAC=Bac_G(3)
                  FV2=Bac_G(3)*N2cBAC(ng)
                  FV3=Bac_G(3)*Fe2cBAC(ng)
                  Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))+    &
     &                                     FV2
                  Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))+    &
     &                                     (NupDOP_ba(i,k,ibac)+        &
     &                                      NupPO4_ba(i,k,ibac))
                  Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(3))*I_Bac_Ceff(ng)
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DON from C2nDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c(ng)*                                &
!!   &                DOC_frac(i,k)*                                    &
     &                (Bio(i,k,iDOMN(ilab))/                            &
     &                 Bio(i,k,iDOMC(ilab)))*Frac_ExBAC_n(ng)
                  FV5=FV2-(NupDON_ba(i,k,ibac)+                         &
     &                     NupNH4_ba(i,k,ibac)-FV4)
!
!  If FV5 is positive then released DON is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDON1=FV4
                    NupNH4_ba(i,k,ibac)=NupNH4_ba(i,k,ibac)+FV5
                  ELSE
                    RelDON1=FV4-FV5
                  END IF
                  NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)-RelDON1
!
!  Release Fe.
!
                  RelFe=NupFe_ba(i,k,ibac)-FV3
                  NupFe_ba(i,k,ibac)=FV3
!
!  Fe limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE
                  Het_BAC=Bac_G(4)
                  FV2=Bac_G(4)*N2cBAC(ng)
                  FV3=Bac_G(4)*P2cBAC(ng)
                  Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))+    &
     &                                     FV2
                  Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))+    &
     &                                     FV3
                  Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))+    &
     &                                     NupFe_ba(i,k,ibac)
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(4))*I_Bac_Ceff(ng)
                  NupDOC_ba(i,k,ibac)=NupDOC_ba(i,k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DON from C2nDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c(ng)*                                &
!!   &                DOC_frac(i,k)*                                    &
     &                Bio(i,k,iDOMN(ilab))/                             &
     &                Bio(i,k,iDOMC(ilab))*Frac_ExBAC_n(ng)
                  FV5=FV2-(NupDON_ba(i,k,ibac)+                         &
     &                     NupNH4_ba(i,k,ibac)-FV4)
!
!  If FV5 is positive then released DON is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDON1=FV4
                    NupNH4_ba(i,k,ibac)=NupNH4_ba(i,k,ibac)+FV5
                  ELSE
                    RelDON1=FV4-FV5
                  END IF
                  NupDON_ba(i,k,ibac)=NupDON_ba(i,k,ibac)-RelDON1
!
!  To get accurate DOP from C2pDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c(ng)*                                &
!!   &               DOC_frac(i,k)*                                     &
     &                Bio(i,k,iDOMP(ilab))/                             &
     &                Bio(i,k,iDOMC(ilab))
                  FV5=FV2-(NupDOP_ba(i,k,ibac)+                         &
     &                     NupPO4_ba(i,k,ibac)-FV4)
!
!  If FV5 is positive then released DOP is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDOP1=FV4
                    NupPO4_ba(i,k,ibac)=NupPO4_ba(i,k,ibac)+FV5
                  ELSE
                    RelDOP1=FV4-FV5
                  END IF
                  NupDOP_ba(i,k,ibac)=NupDOP_ba(i,k,ibac)-RelDOP1
                END IF
!
!  Increment nutrient arrays.
!
                Bio_new(i,k,iBacC(ibac))=Bio_new(i,k,iBacC(ibac))+      &
     &                                   Het_BAC
                FV1=NupDOC_ba(i,k,ibac)-Het_BAC
                Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)+                  &
     &                             FV1
!
!  NOTE: to be strictly accurate we should remove RelDOC1 from DOCNP1,
!       and then add it back, since NupDOC_ba is a net term. This should
!       wash out in the budgeting.
!
                Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))-      &
     &                                   (totDOC_d(i,k)-RelDOC1)
!!   &                                   (totDOC_d(i,k)-RelDOC1)*       &
!!   &                                   DOC_frac(i,k)
!!              Bio_new(i,k,iCDMC(ilab))=Bio_new(i,k,iCDMC(ilab))-      &
!!   &                                   (totDOC_d(i,k)-RelDOC1)*       &
!!   &                                   (1.0_r8-DOC_frac(i,k))
!!
!  This is inclusive of RelDOX1, excretion of DON1 removed above.
!
                Bio_new(i,k,iDOMN(ilab))=Bio_new(i,k,iDOMN(ilab))-      &
     &                                   NupDON_ba(i,k,ibac)
                Bio_new(i,k,iDOMP(ilab))=Bio_new(i,k,iDOMP(ilab))-      &
     &                                   NupDOP_ba(i,k,ibac)
                Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)-                  &
     &                             NupNH4_ba(i,k,ibac)
                Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)-                  &
     &                             NupPO4_ba(i,k,ibac)
                Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)-                  &
     &                             NupFe_ba(i,k,ibac)
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton Losses.
!-----------------------------------------------------------------------
!
          DO iphy=1,Nphy
            DO k=1,N(ng)
              DO i=Istr,Iend
!
!  Excretion.
!
                IF ((C2nALG(i,k,iphy).ge.                               &
     &               C2nALGminABS(iphy,ng)).and.                        &
     &              (C2pALG(i,k,iphy).ge.                               &
     &               C2pALGminABS(iphy,ng)).and.                        &
     &              (HsSiO(iphy,ng).gt.LARGE)) THEN
                  FV1=Bio(i,k,iPhyC(iphy))*ExALG(iphy,ng)
                  Bio_new(i,k,iPhyC(iphy))=Bio_new(i,k,iPhyC(iphy))-    &
     &                                     FV1
!
!  No excretion of CDOC.
!
                  Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))+    &
     &                                     FV1
                ELSE IF ((C2nALG(i,k,iphy).ge.                          &
     &                    C2nALGminABS(iphy,ng)).and.                   &
     &                   (C2pALG(i,k,iphy).ge.                          &
     &                    C2pALGminABS(iphy,ng)).and.                   &
     &                   (C2sALG(i,k,iphy).ge.                          &
     &                    C2SiALGminABS(iphy,ng))) THEN
                  FV1=Bio(i,k,iPhyC(iphy))*ExALG(iphy,ng)
                  Bio_new(i,k,iPhyC(iphy))=Bio_new(i,k,iPhyC(iphy))-    &
     &                                     FV1
!
!  No excretion of CDOC.
!
                  Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))+    &
     &                                     FV1
                END IF
!
!  Grazing.
!
                IF (Bio(i,k,iPhyC(iphy)).gt.refuge(i,k,iphy)) THEN
!
!  Carbon calculations.
!
                  FV1=graz_act(i,k,iphy)*Bio(i,k,iPhyC(iphy))
                  Bio_new(i,k,iPhyC(iphy))=Bio_new(i,k,iPhyC(iphy))-    &
     &                                     FV1
                  Bio_new(i,k,iFecC(isfc))=Bio_new(i,k,iFecC(isfc))+    &
     &                                     FecPEL(iphy,isfc,ng)*FV1
                  Bio_new(i,k,iFecC(iffc))=Bio_new(i,k,iFecC(iffc))+    &
     &                                     FecPEL(iphy,iffc,ng)*FV1
                  FV3=FecDOC(iphy,ng)*FV1
                  Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))+    &
     &                                     (1.0_r8-cDOCfrac_c(ilab,ng))*&
     &                                     FV3
                  Bio_new(i,k,iCDMC(ilab))=Bio_new(i,k,iCDMC(ilab))+    &
     &                                     cDOCfrac_c(ilab,ng)*FV3
                  Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)+                &
     &                               FecCYC(iphy,ng)*FV1
!
!  Nitrogen calculations.
!
                  FV2=graz_act(i,k,iphy)*Bio(i,k,iPhyN(iphy))
                  Bio_new(i,k,iPhyN(iphy))=Bio_new(i,k,iPhyN(iphy))-    &
     &                                     FV2
                  Bio_new(i,k,iFecN(isfc))=Bio_new(i,k,iFecN(isfc))+    &
     &                                     FecPEL(iphy,isfc,ng)*FV2
                  Bio_new(i,k,iFecN(iffc))=Bio_new(i,k,iFecN(iffc))+    &
     &                                     FecPEL(iphy,iffc,ng)*FV2
                  Bio_new(i,k,iDOMN(ilab))=Bio_new(i,k,iDOMN(ilab))+    &
     &                                     FecDOC(iphy,ng)*FV2
                  Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)+                &
     &                               FecCYC(iphy,ng)*FV2
!
!  Silica calculations.
!
                  IF (iPhyS(iphy).gt.0) THEN
                    FV2=graz_act(i,k,iphy)*Bio(i,k,iPhyS(iphy))
                    Bio_new(i,k,iPhyS(iphy))=Bio_new(i,k,iPhyS(iphy))-  &
     &                                       FV2
!
!  Assuming that the fraction of material lost via sloppy feeding/cell
!  lysis also results in silica tests being put into FecS pool.
!
                    Bio_new(i,k,iFecS(isfc))=Bio_new(i,k,iFecS(isfc))+  &
     &                                       FecDOC(iphy,ng)*FV2
                    Bio_new(i,k,iFecS(iffc))=Bio_new(i,k,iFecS(iffc))+  &
     &                                       (1.0_r8-FecDOC(iphy,ng))*  &
     &                                       FV2
                  END IF
!
!  Phosphorus calculations.
!
                  FV2=graz_act(i,k,iphy)*Bio(i,k,iPhyP(iphy))
                  Bio_new(i,k,iPhyP(iphy))=Bio_new(i,k,iPhyP(iphy))-    &
     &                                     FV2
                  Bio_new(i,k,iFecP(isfc))=Bio_new(i,k,iFecP(isfc))+    &
     &                                     FecPEL(iphy,isfc,ng)*FV2
                  Bio_new(i,k,iFecP(iffc))=Bio_new(i,k,iFecP(iffc))+    &
     &                                     FecPEL(iphy,iffc,ng)*FV2
                  Bio_new(i,k,iDOMP(ilab))=Bio_new(i,k,iDOMP(ilab))+    &
     &                                     FecDOC(iphy,ng)*FV2
                  Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)+                &
     &                               FecCYC(iphy,ng)*FV2
!
!  Iron calculations. Assuming no DOMF.
!
                  FV2=graz_act(i,k,iphy)*Bio(i,k,iPhyF(iphy))
                  Bio_new(i,k,iPhyF(iphy))=Bio_new(i,k,iPhyF(iphy))-    &
     &                                     FV2
                  Bio_new(i,k,iFecF(isfc))=Bio_new(i,k,iFecF(isfc))+    &
     &                                     FecPEL(iphy,isfc,ng)*FV2
                  Bio_new(i,k,iFecF(iffc))=Bio_new(i,k,iFecF(iffc))+    &
     &                                     FecPEL(iphy,iffc,ng)*FV2
                  Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)+                &
     &                               (FecCYC(iphy,ng)+                  &
     &                                FecDOC(iphy,ng))*FV2
                END IF
              END DO
            END DO
          END DO
!
!  Pigment Grazing.  No fecal or dissolved terms for pigments.
!
          DO ipig=1,Npig
            DO iphy=1,Nphy
              IF (iPigs(iphy,ipig).gt.0) THEN
                itrc=iPigs(iphy,ipig)
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    IF (Bio(i,k,iPhyC(iphy)).gt.refuge(i,k,iphy)) THEN
                      FV1=graz_act(i,k,iphy)*Bio(i,k,itrc)
                      Bio_new(i,k,itrc)=Bio_new(i,k,itrc) - FV1
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Bacterial losses.
!-----------------------------------------------------------------------
!
!  NOTE: Bacterial growth is completely reminerialized.
!
          DO ibac=1,Nbac
            DO k=1,N(ng)
              DO i=Istr,Iend
!
!  Grazing calculation. (All fecal material to slow sinking pool.)
!
!! WPB - There appears to be some rounding errors that cause bacteria
!!       populations to drop just below initialization values.  Once
!!       they do, they never recover and the new lower values propagate
!!       through the model.  Only evident in Bac_P1 at the moment.
!!
!!                FV1=BacCYC(ng)*Bio_new(i,k,iBacC(ibac))
!!                FV2=BacPEL(ng)*Bio_new(i,k,iBacC(ibac))
!!                FV3=BacDOC(ng)*Bio_new(i,k,iBacC(ibac))
!!                FV4=FV1+FV2+FV3
!
!  Carbon calculations.
!
                Bio_new(i,k,iBacC(ibac))=Bio_new(i,k,iBacC(ibac))-      &
!!   &                                   FV4
     &                                   Bio_new(i,k,iBacC(ibac))
                Bio_new(i,k,iFecC(isfc))=Bio_new(i,k,iFecC(isfc))+      &
!!   &                                   FV2
     &                                   Bio_new(i,k,iBacC(ibac))*      &
     &                                   BacPEL(ng)
                Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))+      &
     &                                   (1.0_r8-cDOCfrac_c(ilab,ng))*  &
!!   &                                   FV3
     &                                   Bio_new(i,k,iBacC(ibac))*      &
     &                                   BacDOC(ng)
                Bio_new(i,k,iCDMC(ilab))=Bio_new(i,k,iCDMC(ilab))+      &
     &                                   cDOCfrac_c(ilab,ng)*           &
!!   &                                   FV3
     &                                   Bio_new(i,k,iBacC(ibac))*      &
     &                                   BacDOC(ng)
                Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)+                  &
!!   &                             FV1
     &                             Bio_new(i,k,iBacC(ibac))*            &
     &                             BacCYC(ng)
!
!  Nitrogen calculations.
!
                Bio_new(i,k,iBacN(ibac))=Bio_new(i,k,iBacN(ibac))-      &
!!   &                                   N2cBAC(ng)*FV4
     &                                   Bio_new(i,k,iBacN(ibac))
                Bio_new(i,k,iFecN(isfc))=Bio_new(i,k,iFecN(isfc))+      &
!!   &                                   N2cBAC(ng)*FV2
     &                                   Bio_new(i,k,iBacN(ibac))*      &
     &                                   BacPEL(ng)
                Bio_new(i,k,iDOMN(ilab))=Bio_new(i,k,iDOMN(ilab))+      &
!!   &                                   N2cBAC(ng)*FV3
     &                                   Bio_new(i,k,iBacN(ibac))*      &
     &                                   BacDOC(ng)
                Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)+                  &
!!   &                             N2cBAC(ng)*FV1
     &                             Bio_new(i,k,iBacN(ibac))*            &
     &                             BacCYC(ng)
!
!  Phosphorous calculations.
!
                Bio_new(i,k,iBacP(ibac))=Bio_new(i,k,iBacP(ibac))-      &
!!   &                                   P2cBAC(ng)*FV4
     &                                   Bio_new(i,k,iBacP(ibac))
                Bio_new(i,k,iFecP(isfc))=Bio_new(i,k,iFecP(isfc))+      &
!!   &                                   P2cBAC(ng)*FV2
     &                                   Bio_new(i,k,iBacP(ibac))*      &
     &                                   BacPEL(ng)
                Bio_new(i,k,iDOMP(ilab))=Bio_new(i,k,iDOMP(ilab))+      &
!!   &                                   P2cBAC(ng)*FV3
     &                                   Bio_new(i,k,iBacP(ibac))*      &
     &                                   BacDOC(ng)
                Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)+                  &
!!   &                             P2cBAC(ng)*FV1
     &                             Bio_new(i,k,iBacP(ibac))*            &
     &                             BacCYC(ng)
!
!  Iron calculations.
!
                Bio_new(i,k,iBacF(ibac))=Bio_new(i,k,iBacF(ibac))-      &
!!   &                                   Fe2cBAC(ng)*FV4
     &                                   Bio_new(i,k,iBacF(ibac))
                Bio_new(i,k,iFecF(isfc))=Bio_new(i,k,iFecF(isfc))+      &
!!   &                                   Fe2cBAC(ng)*FV2
     &                                   Bio_new(i,k,iBacF(ibac))*      &
     &                                   BacPEL(ng)
                Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)+                  &
!!   &                             Fe2cBAC(ng)*(FV1+FV3)
     &                             Bio_new(i,k,iBacF(ibac))*            &
     &                             (BacDOC(ng)+BacCYC(ng))
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Fecal pellet remineralization.
!-----------------------------------------------------------------------
!
          DO ifec=1,Nfec
            DO k=1,N(ng)
              DO i=Istr,Iend
!
!  Carbon calculations.  All carbon goes to CO2.
!
                FV3=Regen_C(i,k,ifec)*Bio(i,k,iFecC(ifec))
                Bio_new(i,k,iFecC(ifec))=Bio_new(i,k,iFecC(ifec))-      &
     &                                   FV3
                Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)+                  &
     &                             FV3
!
!  Nitrogen calculations.  Nitrogen goes to NH4.
!
                FV2=Regen_N(i,k,ifec)*Bio(i,k,iFecN(ifec))
                Bio_new(i,k,iFecN(ifec))=Bio_new(i,k,iFecN(ifec))-      &
     &                                   FV2
                Bio_new(i,k,iNH4_)=Bio_new(i,k,iNH4_)+                  &
     &                             FV2
!
!  Silica calculations.
!
                FV2=Regen_S(i,k,ifec)*Bio(i,k,iFecS(ifec))
                Bio_new(i,k,iFecS(ifec))=Bio_new(i,k,iFecS(ifec))-      &
     &                                   FV2
                Bio_new(i,k,iSiO_)=Bio_new(i,k,iSiO_)+                  &
     &                             FV2
!
!  Phosphorous calculations.
!
                FV2=Regen_P(i,k,ifec)*Bio(i,k,iFecP(ifec))
                Bio_new(i,k,iFecP(ifec))=Bio_new(i,k,iFecP(ifec))-      &
     &                                   FV2
                Bio_new(i,k,iPO4_)=Bio_new(i,k,iPO4_)+                  &
     &                             FV2
!
!  Iron calculations.
!
                FV2=Regen_F(i,k,ifec)*Bio(i,k,iFecF(ifec))
                Bio_new(i,k,iFecF(ifec))=Bio_new(i,k,iFecF(ifec))-      &
     &                                   FV2
                Bio_new(i,k,iFeO_)=Bio_new(i,k,iFeO_)+                  &
     &                             FV2
              END DO
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  CDMC photolysis calculations.
!-----------------------------------------------------------------------
!
          IF (RtUVR_flag(ng)) THEN
            DO i=Istr,Iend
!
!  If Ed_nz(i,N(ng)) > zero, then there is sunlight. Standardizing rate
!  to 1500 umol quanta m-2 s-1.
!
              IF (Ed_nz(i,N(ng)).ge.0.01) THEN
!
                FV1=RtUVR_DIC(ng)*Ed_nz(i,N(ng))/1500.0_r8
                FV2=RtUVR_DOC(ng)*Ed_nz(i,N(ng))/1500.0_r8

!
!  FV4 equals the CDMC1 absorption at 410 nm. 0.012 converts to g m-3.
!  FV5 equals the CDMC2 absorption at 410 nm.
!  Weighted average attenuation of UVB of water at 300 nm = 0.2 m-1.
!
                FV4=Bio(i,N(ng),iCDMC(ilab))*0.012_r8*aDOC410(ilab)
                FV5=Bio(i,N(ng),iCDMC(irct))*0.012_r8*aDOC410(irct)
                photo_decay=0.5_r8*Hz(i,j,N(ng))*                       &
     &                      (0.2_r8+(FV4+FV5)*aDOC300(ilab))
                FV3=EXP(-photo_decay)
                photo_decay=2.0_r8*photo_decay
!
!  Do not photolyze below the euphotic zone.
!
                DO k=N(ng),Keuphotic(i),-1
                  IF (FV3.gt.0.01_r8) THEN
                    FV6=FV5+FV4
                    IF (FV6.gt.0.0_r8) THEN
                      FV7=FV4/FV6
                      photo_DIC=FV3*FV1*FV6
                      photo_DOC=FV3*FV2*FV6
                      total_photo=photo_DIC+photo_DOC
!
!  NOTE: not testing for excess photolysis (CDOC going negative).
!
                      FV4=(1.0_r8-FV7)*total_photo
                      Bio_new(i,k,iCDMC(irct))=Bio_new(i,k,iCDMC(irct))-&
     &                                         FV4
                      Bio_new(i,k,iDOMC(ilab))=Bio_new(i,k,iDOMC(ilab))+&
     &                                         photo_DOC
                      Bio_new(i,k,iCDMC(ilab))=Bio_new(i,k,iCDMC(ilab))-&
     &                                         FV7*total_photo
                      Bio_new(i,k,iDIC_)=Bio_new(i,k,iDIC_)+            &
     &                                         photo_DIC
                    END IF
!
!  FV4 equals the CDMC1 absorption at 410 nm. 0.012 converts to g m-3.
!  FV5 equals the CDMC2 absorption at 410 nm.
!  Weighted average attenuation of UVB of water at 300 nm = 0.2 m-1.
!
                    FV4=Bio(i,k,iCDMC(ilab))*0.012_r8*aDOC410(ilab)
                    FV5=Bio(i,k,iCDMC(irct))*0.012_r8*aDOC410(irct)
                    FV7=photo_decay+                                    &
     &                 0.5_r8*Hz(i,j,k)*(0.2_r8+(FV4+FV5)*aDOC300(ilab))
!
!  If k is greater than the bottom of the euphotic zone (and by
!  by extension the bottom boundary) or the decay constant is
!  greater than 4.61 (or < 1% photolysis zone) then exit do loop.
!
                    FV3=EXP(-FV7)
!
!  Store value for passage through entire Hz(i,j,k).
!
                    photo_decay=photo_decay+2.0_r8*FV7
                  END IF
                END DO
              END IF
            END DO
          END IF
!
!-----------------------------------------------------------------------
!  Create optimal pigment ratios.
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
            IF (Keuphotic(i).le.N(ng)) THEN
              DO iphy=1,Nphy
!
!  Carbon to chlorophyll a ratio
!  This statement says that nutrient limitation of C2CHL ratio overides
!  light adaptation. Minimum of two functions may be more ecologically
!  accurate?
!
                DO k=N(ng),Keuphotic(i),-1
                  IF (b_C2Cn(iphy,ng).lt.0.0_r8+SMALL) THEN
                    C2CHL_w(k,iphy)=MIN((b_C2Cl(iphy,ng)+               &
     &                                   mxC2Cl(iphy,ng)*E0_nz(i,k)),   &
     &                                  C2CHL_max(iphy,ng))
                  ELSE IF (C2nALG(i,k,iphy).gt.                         &
     &                     minC2nALG(iphy,ng)+SMALL) THEN
                    C2CHL_w(k,iphy)=b_C2Cn(iphy,ng)+                    &
     &                              mxC2Cn(iphy,ng)*                    &
     &                              (C2nALG(i,k,iphy)-                  &
     &                               minC2nALG(iphy,ng))
                  ELSE
                    C2CHL_w(k,iphy)=MIN((b_C2Cl(iphy,ng)+               &
     &                                   mxC2Cl(iphy,ng)*E0_nz(i,k)),   &
     &                                  C2CHL_max(iphy,ng))
                  END IF
                END DO
!
!  Chlorophyll a concentation per species. form g CHL a / g C
!
                DO k=N(ng),Keuphotic(i),-1
                  Pigs_w(k,iphy,ichl)=1.0_r8/C2CHL_w(k,iphy)
                END DO
!
!  Chlorophyll b concentration per species. form g CHL b / g C
!
                IF (iPigs(iphy,2).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,2)=b_ChlB(iphy,ng)+                   &
     &                               mxChlB(iphy,ng)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,2)=Pigs_w(k,iphy,2)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Chlorophyll c concentration per species. form g CHL c / g C
!
                IF (iPigs(iphy,3).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,3)=b_ChlC(iphy,ng)+                   &
     &                               mxChlC(iphy,ng)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                 b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,3)=Pigs_w(k,iphy,3)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Photosynthetic caroteniods per species. form g PSC / g C
!
                IF (iPigs(iphy,4).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,4)=b_PSC(iphy,ng)+                    &
     &                               mxPSC(iphy,ng)*                    &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,4)=Pigs_w(k,iphy,4)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Photoprotective caroteniods per species. form g PPC / g C
!
                IF (iPigs(iphy,5).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,5)=b_PPC(iphy,ng)+                    &
     &                               mxPPC(iphy,ng)*                    &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,5)=Pigs_w(k,iphy,5)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Low Urobilin Phycoeurythin concentration per species. g LPUB / g C
!
                IF (iPigs(iphy,6).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,6)=b_LPUb(iphy,ng)+                   &
     &                               mxLPUb(iphy,ng)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,6)=Pigs_w(k,iphy,6)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  High Urobilin Phycoeurythin concentration per species (g HPUB / g C).
!
                IF (iPigs(iphy,7).gt.0) THEN
                  DO k=N(ng),Keuphotic(i),-1
                    Pigs_w(k,iphy,7)=b_HPUb(iphy,ng)+                   &
     &                               mxHPUb(iphy,ng)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy,ng))
                    Pigs_w(k,iphy,7)=Pigs_w(k,iphy,7)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
              END DO
!
!  Calculate pigment ratio changes.
!  NOTE: 12 factor to convert to ugrams (mg m-3)
!
              DO ipig=1,Npig
                DO iphy=1,Nphy
                  IF (iPigs(iphy,ipig).gt.0) THEN
                    itrc=iPigs(iphy,ipig)
                    DO k=N(ng),Keuphotic(i),-1
                      IF ((Bio(i,k,iPhyC(iphy)).gt.0.0_r8).and.         &
     &                    (Bio(i,k,itrc).gt.0.0_r8)) THEN
                        FV1=Bio(i,k,iPhyC(iphy))*12.0_r8
                        FV2=GtALG_r(i,k,iphy)
                        FV3=FV1/                                        &
     &                      (FV2/Pigs_w(k,iphy,ipig)+                   &
     &                       FV1*(1.0_r8-FV2)/                          &
     &                       Bio(i,k,itrc))
                        Bio_new(i,k,itrc)=Bio_new(i,k,itrc)+            &
     &                                    (FV3-Bio(i,k,itrc))
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END IF
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
            itrc=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,itrc)
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
            cff=dtbio*ABS(Wbio(isink))
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
                Bio(i,k,itrc)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
            DO ifec=1,Nfec
              IF (itrc.eq.iFecN(ifec)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iNO3_)=Bio(i,1,iNO3_)+cff1
                END DO
              ELSE IF (itrc.eq.iFecC(ifec)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iDIC_)=Bio(i,1,iDIC_)+cff1
                END DO
              ELSE IF (itrc.eq.iFecP(ifec)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                END DO
              ELSE IF (itrc.eq.iFecS(ifec)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iSiO_)=Bio(i,1,iSiO_)+cff1
                END DO
              ELSE IF (itrc.eq.iFecF(ifec)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iFeO_)=Bio(i,1,iFeO_)+cff1
                END DO
              END IF
            END DO
            DO iphy=1,Nphy
              IF (itrc.eq.iPhyN(iphy)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iNO3_)=Bio(i,1,iNO3_)+cff1
                END DO
              ELSE IF (itrc.eq.iPhyC(iphy)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iDIC_)=Bio(i,1,iDIC_)+cff1
                END DO
              ELSE IF (itrc.eq.iPhyP(iphy)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                END DO
              ELSE IF (itrc.eq.iPhyS(iphy)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iSiO_)=Bio(i,1,iSiO_)+cff1
                END DO
              ELSE IF (itrc.eq.iPhyF(iphy)) THEN
                DO i=Istr,Iend
                  cff1=FC(i,0)*Hz_inv(i,1)
                  Bio(i,1,iFeO_)=Bio(i,1,iFeO_)+cff1
                END DO
              END IF
            END DO
#endif
          END DO SINK_LOOP
!
!-----------------------------------------------------------------------
!  Update the tendency arrays
!-----------------------------------------------------------------------
!
          DO ibio=1,NBT
            itrc=idbio(ibio)
            DO k=1,N(ng)
              DO i=Istr,Iend
                  Bio(i,k,itrc)=Bio(i,k,itrc)+dtbio*Bio_new(i,k,itrc)
              END DO
            END DO
          END DO

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
        DO ibio=1,NBT
          itrc=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
               cff=Bio(i,k,itrc)-Bio_old(i,k,itrc)
               t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
