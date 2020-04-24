      SUBROUTINE tl_biology (ng,tile)
!
!svn $Id: tl_npzd_Franks.h 995 2020-01-10 04:01:28Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2020 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Nutrient-Phytoplankton-Zooplankton-Detritus Model.                  !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Franks et al, 1986: Behavior of simple plankton model with        !
!      food-level acclimation by herbivores, Marine Biology, 91,       !
!      121-129.                                                        !
!                                                                      !
!***********************************************************************
!
      USE mod_param
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
      IF (Lbiofile(iTLM)) THEN
#else
      IF (Lbiofile(iTLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iTLM)=.FALSE.
        BIONAME(iTLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iTLM, 15, __LINE__, __FILE__)
#endif
      CALL tl_biology_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % tl_z_r,                          &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % tl_z_w,                          &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)

#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 15, __LINE__, __FILE__)
#endif
      RETURN
      END SUBROUTINE tl_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE tl_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            Hz, tl_Hz,                            &
     &                            z_r, tl_z_r,                          &
     &                            z_w, tl_z_w,                          &
     &                            t, tl_t)
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
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 1

      integer :: Iter, i, ibio, isink, itrc, itrmx, j, k, ks
      integer :: Iteradj

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-16_r8

      real(r8) :: cff, cff1, cff2, cff3, dtdays
      real(r8) :: tl_cff, tl_cff1
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: tl_cffL, tl_cffR, tl_cu, tl_dltL, tl_dltR

      real(r8), dimension(Nsink) :: Wbio
      real(r8), dimension(Nsink) :: tl_Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio1
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_qc

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
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iSDet                 ! Small detritus
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wDet(ng)                ! Small detritus
      tl_Wbio(1)=tl_wDet(ng)          ! Small detritus
!
      J_LOOP : DO j=Jstr,Jend
!
!  Compute inverse thickness to avoid repeated divisions.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
            tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
            tl_Hz_inv2(i,k)=-Hz_inv2(i,k)*Hz_inv2(i,k)*                 &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
            tl_Hz_inv3(i,k)=-Hz_inv3(i,k)*Hz_inv3(i,k)*                 &
     &                      (tl_Hz(i,j,k-1)+tl_Hz(i,j,k)+               &
     &                       tl_Hz(i,j,k+1))
          END DO
        END DO
!
!  Clear tl_Bio and Bio arrays.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,ibio)=0.0_r8
              Bio1(i,k,ibio)=0.0_r8
              tl_Bio(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping. In this routine, this is not a problem because
!  we only use index nstp in the right-hand-side equations.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=t(i,j,k,nstp,ibio)
              tl_Bio_old(i,k,ibio)=tl_t(i,j,k,nstp,ibio)
            END DO
          END DO
        END DO
!
!  Determine Correction for negativity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff1=MAX(0.0_r8,eps-Bio_old(i,k,iNO3_))+                    &
     &           MAX(0.0_r8,eps-Bio_old(i,k,iPhyt))+                    &
     &           MAX(0.0_r8,eps-Bio_old(i,k,iZoop))+                    &
     &           MAX(0.0_r8,eps-Bio_old(i,k,iSDet))
            tl_cff1=-(0.5_r8-SIGN(0.5_r8,Bio_old(i,k,iNO3_)-eps))*      &
     &               tl_Bio_old(i,k,iNO3_)-                             &
     &               (0.5_r8-SIGN(0.5_r8,Bio_old(i,k,iPhyt)-eps))*      &
     &               tl_Bio_old(i,k,iPhyt)-                             &
     &               (0.5_r8-SIGN(0.5_r8,Bio_old(i,k,iZoop)-eps))*      &
     &               tl_Bio_old(i,k,iZoop)-                             &
     &               (0.5_r8-SIGN(0.5_r8,Bio_old(i,k,iSDet)-eps))*      &
     &               tl_Bio_old(i,k,iSDet)
!
!  If correction needed, determine the largest pool to debit.
!
            IF (cff1.gt.0.0) THEN
              itrmx=idbio(1)
              cff=t(i,j,k,nstp,itrmx)
              DO ibio=idbio(2),idbio(NBT)
                IF (t(i,j,k,nstp,ibio).gt.cff) THEN
                  itrmx=ibio
                  cff=t(i,j,k,nstp,ibio)
                END IF
              END DO
!
!  Update new values.
!
              DO itrc=1,NBT
                ibio=idbio(itrc)
                Bio(i,k,ibio)=MAX(eps,Bio_old(i,k,ibio))-               &
     &                        cff1*                                     &
     &                        (SIGN(0.5_r8, REAL(itrmx-ibio,r8)**2)+    &
     &                         SIGN(0.5_r8,-REAL(itrmx-ibio,r8)**2))
                tl_Bio(i,k,ibio)=(0.5_r8-                               &
     &                            SIGN(0.5_r8,eps-Bio_old(i,k,ibio)))*  &
     &                           tl_Bio_old(i,k,ibio)-                  &
     &                           tl_cff1*                               &
     &                           (SIGN(0.5_r8, REAL(itrmx-ibio,r8)**2)+ &
     &                            SIGN(0.5_r8,-REAL(itrmx-ibio,r8)**2))
              END DO
            ELSE
              DO itrc=1,NBT
                ibio=idbio(itrc)
                Bio(i,k,ibio)=Bio_old(i,k,ibio)
                tl_Bio(i,k,ibio)=tl_Bio_old(i,k,ibio)
              END DO
            END IF
          END DO
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
!  Compute appropriate basic state arrays I.
!
!  Determine Correction for negativity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=MAX(0.0_r8,eps-Bio_old(i,k,iNO3_))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iPhyt))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iZoop))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iSDet))
!
!  If correction needed, determine the largest pool to debit.
!
              IF (cff1.gt.0.0) THEN
                itrmx=idbio(1)
                cff=t(i,j,k,nstp,itrmx)
                DO ibio=idbio(2),idbio(NBT)
                  IF (t(i,j,k,nstp,ibio).gt.cff) THEN
                    itrmx=ibio
                    cff=t(i,j,k,nstp,ibio)
                  END IF
                END DO
!
!  Update new values.
!
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=MAX(eps,Bio_old(i,k,ibio))-             &
     &                          cff1*(SIGN(0.5_r8,                      &
     &                                      REAL(itrmx-ibio,r8)**2)+    &
     &                                SIGN(0.5_r8,                      &
     &                                     -REAL(itrmx-ibio,r8)**2))
                END DO
              ELSE
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=Bio_old(i,k,ibio)
                END DO
              END IF
            END DO
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Nutrient uptake by phytoplankton.
!
            cff1=dtdays*Vm_NO3(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iPhyt)*                                     &
     &              cff1*EXP(K_ext(ng)*z_r(i,j,k))/                     &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
                Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/                          &
     &                         (1.0_r8+cff)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!  Phytoplankton grazing by Zooplankton and mortality to Detritus
!  (rate: PhyMR).
!
              cff1=dtdays*ZooGR(ng)
              cff2=dtdays*PhyMR(ng)
              cff3=K_phy(ng)*K_phy(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=Bio(i,k,iZoop)*Bio(i,k,iPhyt)*cff1/               &
     &                (cff3+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                  Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/                        &
     &                           (1.0_r8+cff+cff2)
                  Bio1(i,k,iZoop)=Bio(i,k,iZoop)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)+                        &
     &                           Bio(i,k,iPhyt)*cff*(1.0_r8-ZooGA(ng))
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)+                        &
     &                           Bio(i,k,iPhyt)*                        &
     &                           (cff2+cff*(ZooGA(ng)-ZooEC(ng)))
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iPhyt)*cff*ZooEC(ng)
                END DO
              END DO
!
!  Zooplankton excretion to nutrients and mortality to Detritus.
!
              cff1=1.0_r8/(1.0_r8+dtdays*(ZooMR(ng)+ZooMD(ng)))
              cff2=dtdays*ZooMR(ng)
              cff3=dtdays*ZooMD(ng)
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
!  Detritus breakdown to nutrients.
!
              cff1=dtdays*DetRR(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff2
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iSDet)*cff1
                END DO
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
!  Nutrient uptake by phytoplankton.
!
          cff1=dtdays*Vm_NO3(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio1(i,k,iPhyt)*                                      &
     &            cff1*EXP(K_ext(ng)*z_r(i,j,k))/                       &
     &            (K_NO3(ng)+Bio1(i,k,iNO3_))
              tl_cff=(tl_Bio(i,k,iPhyt)*                                &
     &                cff1*EXP(K_ext(ng)*z_r(i,j,k))-                   &
     &                tl_Bio(i,k,iNO3_)*cff)/                           &
     &               (K_NO3(ng)+Bio1(i,k,iNO3_))+                       &
     &               K_ext(ng)*tl_z_r(i,j,k)*cff
!>            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/                            &
!>   &                       (1.0_r8+cff)
!>
              tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-                     &
     &                           tl_cff*Bio(i,k,iNO3_))/                &
     &                          (1.0_r8+cff)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                            &
!>   &                       Bio(i,k,iNO3_)*cff
!>
              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+                      &
     &                          tl_Bio(i,k,iNO3_)*cff+                  &
     &                          Bio(i,k,iNO3_)*tl_cff
            END DO
          END DO
!
!  Compute appropriate basic state arrays II.
!
!  Determine Correction for negativity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=MAX(0.0_r8,eps-Bio_old(i,k,iNO3_))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iPhyt))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iZoop))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iSDet))
!
!  If correction needed, determine the largest pool to debit.
!
              IF (cff1.gt.0.0) THEN
                itrmx=idbio(1)
                cff=t(i,j,k,nstp,itrmx)
                DO ibio=idbio(2),idbio(NBT)
                  IF (t(i,j,k,nstp,ibio).gt.cff) THEN
                    itrmx=ibio
                    cff=t(i,j,k,nstp,ibio)
                  END IF
                END DO
!
!  Update new values.
!
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=MAX(eps,Bio_old(i,k,ibio))-             &
     &                          cff1*(SIGN(0.5_r8,                      &
     &                                      REAL(itrmx-ibio,r8)**2)+    &
     &                                SIGN(0.5_r8,                      &
     &                                     -REAL(itrmx-ibio,r8)**2))
                END DO
              ELSE
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=Bio_old(i,k,ibio)
                END DO
              END IF
            END DO
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Nutrient uptake by phytoplankton.
!
            cff1=dtdays*Vm_NO3(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iPhyt)*                                     &
     &              cff1*EXP(K_ext(ng)*z_r(i,j,k))/                     &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
                Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/                          &
     &                         (1.0_r8+cff)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff
              END DO
            END DO
!
!  Phytoplankton grazing by Zooplankton and mortality to Detritus
!  (rate: PhyMR).
!
            cff1=dtdays*ZooGR(ng)
            cff2=dtdays*PhyMR(ng)
            cff3=K_phy(ng)*K_phy(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iZoop)*Bio(i,k,iPhyt)*cff1/                 &
     &              (cff3+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/                          &
     &                         (1.0_r8+cff+cff2)
                Bio1(i,k,iZoop)=Bio(i,k,iZoop)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         Bio(i,k,iPhyt)*cff*(1.0_r8-ZooGA(ng))
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*                          &
     &                         (cff2+cff*(ZooGA(ng)-ZooEC(ng)))
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*cff*ZooEC(ng)
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!  Zooplankton excretion to nutrients and mortality to Detritus.
!
              cff1=1.0_r8/(1.0_r8+dtdays*(ZooMR(ng)+ZooMD(ng)))
              cff2=dtdays*ZooMR(ng)
              cff3=dtdays*ZooMD(ng)
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
!  Detritus breakdown to nutrients.
!
              cff1=dtdays*DetRR(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff2
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                        &
     &                           Bio(i,k,iSDet)*cff1
                END DO
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
!  Phytoplankton grazing by Zooplankton and mortality to Detritus
!  (rate: PhyMR).
!
          cff1=dtdays*ZooGR(ng)
          cff2=dtdays*PhyMR(ng)
          cff3=K_phy(ng)*K_phy(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio1(i,k,iZoop)*Bio1(i,k,iPhyt)*cff1/                 &
     &            (cff3+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              tl_cff=((tl_Bio(i,k,iZoop)*Bio1(i,k,iPhyt)+               &
     &                 Bio1(i,k,iZoop)*tl_Bio(i,k,iPhyt))*cff1-         &
     &                2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)*cff)/    &
     &               (cff3+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/                            &
!>   &                       (1.0_r8+cff+cff2)
!>
              tl_Bio(i,k,iPhyt)=(tl_Bio(i,k,iPhyt)-                     &
     &                           tl_cff*Bio(i,k,iPhyt))/                &
     &                          (1.0_r8+cff+cff2)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)+                            &
!>   &                       Bio(i,k,iPhyt)*cff*(1.0_r8-ZooGA(ng))
!>
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)+                      &
     &                          tl_Bio(i,k,iPhyt)*                      &
     &                          cff*(1.0_r8-ZooGA(ng))+                 &
     &                          Bio(i,k,iPhyt)*                         &
     &                          tl_cff*(1.0_r8-ZooGA(ng))
!>            Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
!>   &                       Bio(i,k,iPhyt)*                            &
!>   &                       (cff2+cff*(ZooGA(ng)-ZooEC(ng)))
!>
              tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)+                      &
     &                          tl_Bio(i,k,iPhyt)*                      &
     &                          (cff2+cff*(ZooGA(ng)-ZooEC(ng)))+       &
     &                          Bio(i,k,iPhyt)*                         &
     &                          tl_cff*(ZooGA(ng)-ZooEC(ng))
!>            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
!>   &                       Bio(i,k,iPhyt)*cff*ZooEC(ng)
!>
              tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
     &                          tl_Bio(i,k,iPhyt)*cff*ZooEC(ng)+        &
     &                          Bio(i,k,iPhyt)*tl_cff*ZooEC(ng)
            END DO
          END DO
!
!  Zooplankton excretion to nutrients and mortality to Detritus.
!
          cff1=1.0_r8/(1.0_r8+dtdays*(ZooMR(ng)+ZooMD(ng)))
          cff2=dtdays*ZooMR(ng)
          cff3=dtdays*ZooMD(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
!>
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)*cff1
!>            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
!>   &                       Bio(i,k,iZoop)*cff2
!>
              tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
     &                          tl_Bio(i,k,iZoop)*cff2
!>            Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
!>   &                       Bio(i,k,iZoop)*cff3
!>
              tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)+                      &
     &                          tl_Bio(i,k,iZoop)*cff3
            END DO
          END DO
!
!  Detritus breakdown to nutrients.
!
          cff1=dtdays*DetRR(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff2
!>
              tl_Bio(i,k,iSDet)=tl_Bio(i,k,iSDet)*cff2
!>            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
!>   &                       Bio(i,k,iSDet)*cff1
!>
              tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+                      &
     &                          tl_Bio(i,k,iSDet)*cff1
            END DO
          END DO
!
!  Compute appropriate basic state arrays III.
!
!  Determine Correction for negativity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=MAX(0.0_r8,eps-Bio_old(i,k,iNO3_))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iPhyt))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iZoop))+                  &
     &             MAX(0.0_r8,eps-Bio_old(i,k,iSDet))
!
!  If correction needed, determine the largest pool to debit.
!
              IF (cff1.gt.0.0) THEN
                itrmx=idbio(1)
                cff=t(i,j,k,nstp,itrmx)
                DO ibio=idbio(2),idbio(NBT)
                  IF (t(i,j,k,nstp,ibio).gt.cff) THEN
                    itrmx=ibio
                    cff=t(i,j,k,nstp,ibio)
                  END IF
                END DO
!
!  Update new values.
!
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=MAX(eps,Bio_old(i,k,ibio))-             &
     &                          cff1*(SIGN(0.5_r8,                      &
     &                                      REAL(itrmx-ibio,r8)**2)+    &
     &                                SIGN(0.5_r8,                      &
     &                                     -REAL(itrmx-ibio,r8)**2))
                END DO
              ELSE
                DO itrc=1,NBT
                  ibio=idbio(itrc)
                  Bio(i,k,ibio)=Bio_old(i,k,ibio)
                END DO
              END IF
            END DO
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Nutrient uptake by phytoplankton.
!
            cff1=dtdays*Vm_NO3(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iPhyt)*                                     &
     &              cff1*EXP(K_ext(ng)*z_r(i,j,k))/                     &
     &              (K_NO3(ng)+Bio(i,k,iNO3_))
                Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/                          &
     &                         (1.0_r8+cff)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         Bio(i,k,iNO3_)*cff
              END DO
            END DO
!
!  Phytoplankton grazing by Zooplankton and mortality to Detritus
!  (rate: PhyMR).
!
            cff1=dtdays*ZooGR(ng)
            cff2=dtdays*PhyMR(ng)
            cff3=K_phy(ng)*K_phy(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=Bio(i,k,iZoop)*Bio(i,k,iPhyt)*cff1/                 &
     &              (cff3+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/                          &
     &                         (1.0_r8+cff+cff2)
                Bio1(i,k,iZoop)=Bio(i,k,iZoop)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         Bio(i,k,iPhyt)*cff*(1.0_r8-ZooGA(ng))
                Bio(i,k,iSDet)=Bio(i,k,iSDet)+                          &
     &                         Bio(i,k,iPhyt)*                          &
     &                         (cff2+cff*(ZooGA(ng)-ZooEC(ng)))
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iPhyt)*cff*ZooEC(ng)
              END DO
            END DO
!
!  Zooplankton excretion to nutrients and mortality to Detritus.
!
            cff1=1.0_r8/(1.0_r8+dtdays*(ZooMR(ng)+ZooMD(ng)))
            cff2=dtdays*ZooMR(ng)
            cff3=dtdays*ZooMD(ng)
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
!  Detritus breakdown to nutrients.
!
            cff1=dtdays*DetRR(ng)
            cff2=1.0_r8/(1.0_r8+cff1)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff2
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                          &
     &                         Bio(i,k,iSDet)*cff1
              END DO
            END DO
!
            IF (Iteradj.ne.Iter) THEN
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
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
                qc(i,k)=Bio(i,k,ibio)
                tl_qc(i,k)=tl_Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
                tl_FC(i,k)=(tl_qc(i,k+1)-tl_qc(i,k))*Hz_inv2(i,k)+      &
     &                     (qc(i,k+1)-qc(i,k))*tl_Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                tl_dltR=tl_Hz(i,j,k)*FC(i,k)+Hz(i,j,k)*tl_FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                tl_dltL=tl_Hz(i,j,k)*FC(i,k-1)+Hz(i,j,k)*tl_FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                tl_cff=tl_Hz(i,j,k-1)+2.0_r8*tl_Hz(i,j,k)+tl_Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                tl_cffR=tl_cff*FC(i,k)+cff*tl_FC(i,k)
                cffL=cff*FC(i,k-1)
                tl_cffL=tl_cff*FC(i,k-1)+cff*tl_FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  tl_dltR=0.0_r8
                  dltL=0.0_r8
                  tl_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                  tl_dltR=tl_cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                  tl_dltL=tl_cffR
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
                tl_cff=(tl_dltR-tl_dltL)*Hz_inv3(i,k)+                  &
     &                 (dltR-dltL)*tl_Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                tl_dltR=tl_dltR-tl_cff*Hz(i,j,k+1)-cff*tl_Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                tl_dltL=tl_dltL+tl_cff*Hz(i,j,k-1)+cff*tl_Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                tl_bR(i,k)=tl_qc(i,k)+tl_dltR
                bL(i,k)=qc(i,k)-dltL
                tl_bL(i,k)=tl_qc(i,k)-tl_dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                tl_WR(i,k)=2.0_r8*(2.0_r8*dltR-dltL)*                   &
     &                            (2.0_r8*tl_dltR-tl_dltL)
                WL(i,k)=(dltR-2.0_r8*dltL)**2
                tl_WL(i,k)=2.0_r8*(dltR-2.0_r8*dltL)*                   &
     &                            (tl_dltR-2.0_r8*tl_dltL)
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                tl_dltL=(0.5_r8-SIGN(0.5_r8,cff-WL(i,k  )))*            &
     &                  tl_WL(i,k  )
                dltR=MAX(cff,WR(i,k+1))
                tl_dltR=(0.5_r8-SIGN(0.5_r8,cff-WR(i,k+1)))*            &
     &                  tl_WR(i,k+1)
                bR1(i,k)=bR(i,k)
                bL1(i,k+1)=bL(i,k+1)
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                tl_bR(i,k)=(tl_dltR*bR1(i,k  )+dltR*tl_bR(i,k  )+       &
     &                      tl_dltL*bL1(i,k+1)+dltL*tl_bL(i,k+1))/      &
     &                      (dltR+dltL)-                                &
     &                      (tl_dltR+tl_dltL)*bR(i,k)/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
                tl_bL(i,k+1)=tl_bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
              tl_FC(i,N(ng))=0.0_r8         ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
              tl_bR(i,N(ng))=2.0_r8*tl_qc(i,N(ng))-tl_bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
              tl_bR(i,N(ng))=1.5_r8*tl_qc(i,N(ng))-0.5_r8*tl_bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
              tl_bR(i,N(ng))=tl_qc(i,N(ng)) ! default strictly monotonic
              tl_bL(i,N(ng))=tl_qc(i,N(ng)) ! conditions
              tl_bR(i,N(ng)-1)=tl_qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              tl_bR(i,1)=tl_bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
              tl_bL(i,1)=2.0_r8*tl_qc(i,1)-tl_bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              tl_bR(i,1)=tl_bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
              tl_bL(i,1)=1.5_r8*tl_qc(i,1)-0.5_r8*tl_bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
              tl_bL(i,2)=tl_qc(i,1)         ! bottom grid boxes are
              tl_bR(i,1)=tl_qc(i,1)         ! re-assumed to be
              tl_bL(i,1)=tl_qc(i,1)         ! piecewise constant.
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
                tl_dltR=tl_bR(i,k)-tl_qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                tl_dltL=tl_qc(i,k)-tl_bL(i,k)
                cffR=2.0_r8*dltR
                tl_cffR=2.0_r8*tl_dltR
                cffL=2.0_r8*dltL
                tl_cffL=2.0_r8*tl_dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  tl_dltR=0.0_r8
                  dltL=0.0_r8
                  tl_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                  tl_dltR=tl_cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                  tl_dltL=tl_cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                tl_bR(i,k)=tl_qc(i,k)+tl_dltR
                bL(i,k)=qc(i,k)-dltL
                tl_bL(i,k)=tl_qc(i,k)-tl_dltL
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
            tl_cff=dtdays*SIGN(1.0_r8,Wbio(isink))*tl_Wbio(isink)
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                tl_FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                tl_WL(i,k)=tl_z_w(i,j,k-1)+tl_cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                tl_WR(i,k)=tl_Hz(i,j,k)*qc(i,k)+Hz(i,j,k)*tl_qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                    tl_FC(i,k-1)=tl_FC(i,k-1)+tl_WR(i,ks)
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
                tl_cu=(0.5_r8+SIGN(0.5_r8,                              &
     &                             (1.0_r8-(WL(i,k)-z_w(i,j,ks-1))*     &
     &                             Hz_inv(i,ks))))*                     &
     &                ((tl_WL(i,k)-tl_z_w(i,j,ks-1))*Hz_inv(i,ks)+      &
     &                 (WL(i,k)-z_w(i,j,ks-1))*tl_Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
                tl_FC(i,k-1)=tl_FC(i,k-1)+                              &
     &                       (tl_Hz(i,j,ks)*cu+Hz(i,j,ks)*tl_cu)*       &
     &                       (bL(i,ks)+                                 &
     &                        cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-           &
     &                            (1.5_r8-cu)*                          &
     &                            (bR(i,ks)+bL(i,ks)-                   &
     &                             2.0_r8*qc(i,ks))))+                  &
     &                       Hz(i,j,ks)*cu*                             &
     &                       (tl_bL(i,ks)+                              &
     &                        tl_cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-        &
     &                               (1.5_r8-cu)*                       &
     &                               (bR(i,ks)+bL(i,ks)-                &
     &                               2.0_r8*qc(i,ks)))+                 &
     &                        cu*(0.5_r8*(tl_bR(i,ks)-tl_bL(i,ks))+     &
     &                            tl_cu*                                &
     &                            (bR(i,ks)+bL(i,ks)-2.0_r8*qc(i,ks))-  &
     &                            (1.5_r8-cu)*                          &
     &                            (tl_bR(i,ks)+tl_bL(i,ks)-             &
     &                             2.0_r8*tl_qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                tl_Bio(i,k,ibio)=tl_qc(i,k)+                            &
     &                           (tl_FC(i,k)-tl_FC(i,k-1))*Hz_inv(i,k)+ &
     &                           (FC(i,k)-FC(i,k-1))*tl_Hz_inv(i,k)
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
              tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
!>            t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
!>
              tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
     &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)
            END DO
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE tl_biology_tile
