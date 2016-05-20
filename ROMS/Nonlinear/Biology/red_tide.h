#undef  NEWMAN
#define LIMIT_INTERIOR

      SUBROUTINE biology (ng,tile)
!
!svn $Id: red_tide.h 791 2016-05-05 22:39:42Z arango $
!******************************************************** Ruoying He ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  Red Tide Biological Model: Alexandrium fundyense                    !
!                                                                      !
!  This routine computes the biological sources and sinks and          !
!  updates the global tracer array for dinoflagellates                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Stock, C.A., D.J. McGillicudy, A.R. Solow, and D.M. Anderson,     !
!      2005: Evaluating hypotheses for the initiation and development  !
!      of Alexandrium fundyense blooms in the western Gulf of Maine    !
!      using a coupled physical-biological model, Deep-Sea Research II,!
!      52, 2715-2744.                                                  !
!                                                                      !
!    He, R., D.J. McGillicuddy, B.A. Keafer, and D.M. Anderson, 2008:  !
!      Historic 2005 toxic bloom of Alexandrium fundyense in the       !
!      western Gulf of Maine: 2, Coupled biophysical modeling, J.      !
!      Geophys. Res., 113, C07040, doi:10.1029/2007JC004602.           !
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
#ifdef DAILY_SHORTWAVE
     &                   FORCES(ng) % srflx_avg,                        &
#endif
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % CystIni,                           &
     &                   OCEAN(ng) % DIN_obs,                           &
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
     &                         Hz, z_r, z_w,                            &
#ifdef DAILY_SHORTWAVE
     &                         srflx_avg,                               &
#endif
     &                         srflx, CystIni, DIN_obs,                 &
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
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: CystIni(LBi:,LBj:)
      real(r8), intent(in) :: DIN_obs(LBi:,LBj:,:)
# ifdef DAILY_SHORTWAVE
      real(r8), intent(in) :: srflx_avg(LBi:,LBj:)
# endif
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: CystIni(LBi:UBi,LBj:UBj)
# ifdef DAILY_SHORTWAVE
      real(r8), intent(in) :: srflx_avg(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DIN_ob(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nswim = 1

      integer, parameter :: NsedLayers = 10

      integer :: Iter, i, ibio, iswim, itrc, j, k, ks, ksed
      integer :: iday, month, year

      integer, dimension(Nswim) :: idswim

      real(r8) :: Cell_Flux, C_depth, DIN, E_flux, EndoScale
      real(r8) :: Rad, RadScale
      real(r8) :: GermD, GermL, G_DIN, G_light, G_rate, M_rate
      real(r8) :: G_fac, S_fac, T_fac
      real(r8) :: dtdays, oNsedLayers, salt, temp, wmig
      real(r8) :: hour, yday

      real(r8) :: alpha, cff, cffL, cffR, deltaL, deltaR, dz, wdt

      real(r8), parameter :: eps = 1.0E-8_r8

      real(r8), dimension(Nswim) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: Germ

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: r

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Light
!
!  Temperature dependent grow factor polynomial coefficients.
!
!!    real(r8), parameter :: TC0 = 0.382_r8     ! Stock el al, 2005
!!    real(r8), parameter :: TC1 =-0.0867_r8    ! Eq. 5a
!!    real(r8), parameter :: TC2 = 0.0160_r8
!!    real(r8), parameter :: TC3 =-0.000513_r8
!!
      real(r8), parameter :: TC0 = 0.379_r8     ! Revised
      real(r8), parameter :: TC1 =-0.0961_r8    ! Stock 8/15/2006
      real(r8), parameter :: TC2 = 0.0169_r8
      real(r8), parameter :: TC3 =-0.000536_r8
!
!  Salinity dependent grow factor polynomial coefficients.
!
!!    real(r8), parameter :: SC0 =-0.872_r8     ! Stock el al, 2005
!!    real(r8), parameter :: SC1 = 0.220_r8     ! Eq. 6
!!    real(r8), parameter :: SC2 =-0.00808_r8
!!    real(r8), parameter :: SC3 = 0.0000882_r8
!!
      real(r8), parameter :: SC0 =-0.693_r8     ! Revised
      real(r8), parameter :: SC1 = 0.186_r8     ! Stock 8/15/2006
      real(r8), parameter :: SC2 =-0.00622_r8
      real(r8), parameter :: SC3 = 0.0000557_r8

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
!  Set shortwave radiation scale. In ROMS all the fluxes are kinematic.
!
      RadScale=rho0*Cp          ! Celsius m/s  to  Watts/m2
!
!  Set critical depth (m; negative) used in the growth function.
!
      C_depth=(LOG(G_r(ng)/(G_eff(ng)*srad_Cdepth(ng))))/AttW(ng)
!
!  Set vertical swimming identification vector.
!
      idswim(1)=iDino           ! Dinoflagellate
!
!  Set vertical swimming velocity vector in the same order as the
!  identification vector, IDSWIM.
!
      Wbio(1)=wDino(ng)         ! Dinoflagellate
!
!  Set scale for germination term.
!
      oNsedLayers=1.0_r8/REAL(NsedLayers,r8)
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
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
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),36.0_r8)
            Bio(i,k,isalt)=MAX(0.0_r8,t(i,j,k,nstp,isalt))
          END DO
        END DO
!
!  Calculate endogenous clock scale used in the cysts germination
!  term. The cysts germination rates are regulated by an endogenous
!  circannual clock.
!
        CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
!
        IF (yday.lt.Month_MidDay(1)) THEN
          cff=(365.0_r8-Month_MidDay(12)+yday)/                         &
     &        (365.0_r8-Month_MidDay(12)+Month_MidDay(1))
          EndoScale=GPN(12)+cff*(GPN(1)-GPN(12))
        ELSE IF (yday.ge.Month_MidDay(12)) THEN
          cff=(yday-Month_MidDay(12))/                                  &
              (365.0_r8-Month_MidDay(12)+Month_MidDay(1))
          EndoScale=GPN(12)+cff*(GPN(1)-GPN(12))
        ELSE
          DO i=1,11
            IF ((yday.ge.Month_MidDay(i)).and.                          &
     &          (yday.lt.Month_MidDay(i+1))) THEN
              cff=(yday-Month_MidDay(i))/                               &
     &            (Month_MidDay(i+1)-Month_MidDay(i))
              EndoScale=GPN(i)+cff*(GPN(i+1)-GPN(i))
            END IF
          END DO
        END IF
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Add Cyst germination flux at the bottom layer
!-----------------------------------------------------------------------
!
!  Calculate Cyst germination rate at the top of the sediment layer
!  as a function of bottom water temperature and non-spectral
!  irradiance.
!
          DO i=Istr,Iend
!
!  Calculate "light" and "dark" cyst germination rates as a function
!  of bottom temperature.
!
            temp=Bio(i,1,itemp)                         ! bottom, k=1
            GermL=(1.50_r8+                                              &
     &             (8.72_r8-1.50_r8)*0.5_r8*                             &
     &             (TANH(0.790_r8*temp-6.27_r8)+1.0_r8))*oNsedLayers
            GermD=(1.04_r8+                                              &
     &             (4.26_r8-1.04_r8)*0.5_r8*                             &
     &             (TANH(0.394_r8*temp-3.33_r8)+1.0_r8))*oNsedLayers
!
!  Compute non-spectral irradiance flux at each sediment layer.  Then,
!  compute cyst germination rate according to the light regime.
!
            Germ(i)=0.0_r8                             ! initialize
            DO ksed=1,NsedLayers
# ifdef DAILY_SHORTWAVE
              E_flux=RadScale*srflx_avg(i,j)*                           &
     &               EXP( AttW(ng)*z_w(i,j,0)-                          &
     &                    AttS(ng)*Dg(ng)*(REAL(ksed,r8)-0.5) )
# else
              E_flux=RadScale*srflx(i,j)*                               &
     &               EXP( AttW(ng)*z_w(i,j,0)-                          &
     &                    AttS(ng)*Dg(ng)*(REAL(ksed,r8)-0.5) )
# endif
              IF (E_flux.gt.E_light(ng)) THEN
                Germ(i)=Germ(i)+GermL
              ELSE IF (E_flux.lt.E_dark(ng)) THEN
                Germ(i)=Germ(i)+GermD
              ELSE
                Germ(i)=Germ(i)+                                        &
     &                  (GermD+                                         &
     &                   (GermL-GermD)*                                 &
     &                   ((E_flux-E_dark(ng))/                          &
     &                    (E_light(ng)-E_dark(ng))))
              END IF
            END DO
!
!  Multiply by endogenous clock factor. The cyst germination are
!  regulated by an endogenous circannual clock.
!
            Germ(i)=Germ(i)*EndoScale
!
!  Convert percentage cysts/day into decimal fraction of cysts.
!
            Germ(i)=Germ(i)*0.01_r8
!
!  Calculate the flux of cells away from the bottom. It is referenced
!  to the initial number of cysts to be consistent with laboratory
!  experiments.
!
            Cell_Flux=CystIni(i,j)*                                     &
     &                Germ(i)*Hz_inv(i,1)                     ! cells/m3
!
!  Add cell flux at the bottom layer (k=1).
!
            Bio(i,1,iDino)=Bio(i,1,iDino)+Cell_Flux*dtdays
          END DO
!
!-----------------------------------------------------------------------
!  Compute growth term.
!-----------------------------------------------------------------------
!
!  The growth is dependent on temperature, salinity, non-spectral
!  irradiance (light), and nutrient (Dissolved Inorganic Nutrient, DIN).
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              temp=Bio(i,k,itemp)
              salt=Bio(i,k,isalt)
!
!  Compute Alexandrium fundyense temperature dependent growth factor
!  using a cubic polynomial fitted to available data.
!
              IF (temp.ge.Tmin_growth(ng)) THEN
                T_fac=TC0+temp*(TC1+temp*(TC2+temp*TC3))
              ELSE                               ! linear extrapolation
!!              T_fac=TC0+temp*(TC1+temp*(TC22+temp*TC3))-              &
!!   &                0.0343_r8*(5.0_r8-temp)    ! Stock el al, 2005
!!                                               ! Eq. 5b
                T_fac=0.254_r8-0.0327_r8*(5.0_r8-temp) ! Stock 8/15/2006
              END IF
!
!  Compute Alexandrium fundyense salinity dependent growth factor
!  using a cubic polynomial to fit to available data.
!
              S_fac=SC0+salt*(SC1+salt*(SC2+salt*SC3))
!
!  Compute temperature and salinity growth factor.
!
              G_fac=T_fac*S_fac
!
!  Compute light dependency factor (Platt and Jassby, 1976).
!
# ifdef DAILY_SHORTWAVE
              Rad=srflx_avg(i,j)*RadScale*EXP(AttW(ng)*z_r(i,j,k))
# else
              Rad=srflx(i,j)*RadScale*EXP(AttW(ng)*z_r(i,j,k))
# endif
              IF (z_r(i,j,k).gt.C_depth) THEN
                cff=Gmax(ng)*G_fac+G_r(ng)
                G_light=MAX(0.0_r8,cff*TANH(G_eff(ng)*Rad/cff)-G_r(ng))
              ELSE
                G_light=0.0_r8
              END IF
!
!  Compute dissolved inorganic nutrient (DIN) dependency.
!  [JWilkin: This ELSE block below appears redundant because if
!   z_r(i,j,k).le.C_depth then G_light=0.0_r8 (see above) and
!   therefore G_rate will be set to zero (see below) regardless of
!   the calculated value of G_DIN].
!
              IF (z_r(i,j,k).gt.C_depth) THEN
                DIN=DIN_OBS(i,j,k)
              ELSE
                DIN=DIN_Cdepth(ng)
              END IF
!
!  The nutrient dependence is modeled by the Monod formulation with
!  half-saturation Kn.
!
              G_DIN=Gmax(ng)*G_fac*DIN/(MAX(Kn(ng),0.0_r8)+DIN)
!
!  Compute growth term (implicit).  The growth rate is either limited
!  by the nutrient or light. The rate is capped to be positive.
!
              G_rate=MAX(MIN(G_light,G_DIN),0.0_r8)
              Bio(i,k,iDino)=Bio(i,k,iDino)/(1.0_r8-G_rate*dtdays)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Compute mortality term.
!-----------------------------------------------------------------------
!
!  The mortality is modeled as function dependent on temperature
!  (implicit).  The simple input mortality rate is not used.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              temp=Bio(i,k,itemp)
!!            M_rate=Mor(ng)
              M_rate=0.019_r8+                                          &
     &               0.066_r8*21.76_r8**((temp-10.35_r8)*0.1_r8)
              Bio(i,k,iDino)=Bio(i,k,iDino)/(1.0_r8+M_rate*dtdays)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking/ascending term: dinoflagellate swimming
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,iswim)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to vertical motion.
!  Many thanks to Sasha Shchepetkin for the updated algorithm.
!
          SWIM_LOOP: DO iswim=1,Nswim
            ibio=idswim(iswim)
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(Bio(i,k+1,ibio)-Bio(i,k,ibio))/                &
     &                  (Hz(i,j,k+1)+Hz(i,j,k))
              END DO
            END DO
!
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                deltaR=Hz(i,j,k)*FC(i,k  )
                deltaL=Hz(i,j,k)*FC(i,k-1)
                IF (deltaR*deltaL.lt.0.0_r8) THEN
                  deltaR=0.0_r8
                  deltaL=0.0_r8
                END IF
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k  )
                cffL=cff*FC(i,k-1)
                IF (ABS(deltaR).gt.ABS(cffL)) deltaR=cffL
                IF (ABS(deltaL).gt.ABS(cffR)) deltaL=cffR
                cff=(deltaR-deltaL)/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
                deltaR=deltaR-cff*Hz(i,j,k+1)
                deltaL=deltaL+cff*Hz(i,j,k-1)
!
                aR(i,k)=Bio(i,k,ibio)+deltaR
                aL(i,k)=Bio(i,k,ibio)-deltaL
!
                dR(i,k)=(2.0_r8*deltaR-deltaL)**2
                dL(i,k)=(2.0_r8*deltaL-deltaR)**2
              END DO
            END DO

#ifdef LIMIT_INTERIOR
!
! Apply boundary conditions for strictly monotonic option. The only way
! to avoid extrapolation toward the boundary is to assume that field is
! simply constant within topmost and bottommost grid boxes.
!

            DO i=Istr,Iend
              aR(i,N(ng))=Bio(i,N(ng),ibio)
              aL(i,N(ng))=Bio(i,N(ng),ibio)
              dR(i,N(ng))=0.0_r8
              dL(i,N(ng))=0.0_r8
!
              aR(i,1)=Bio(i,1,ibio)
              aL(i,1)=Bio(i,1,ibio)
              dR(i,1)=0.0_r8
              dL(i,1)=0.0_r8
            END DO
#else
!
! Apply Neumann or linear continuation boundary conditions. Notice that
! for Neumann conditions, the extrapolate values aR(i,N(ng)) and aL(i,0)
! exceed corresponding box values.
!
            DO i=Istr,Iend
              aL(i,N(ng))=aR(i,N(ng)-1)
# ifdef NEUMANN
              aR(i,N(ng))=1.5_r8*Bio(i,N(ng),ibio)-0.5_r8*aL(i,N(ng))
# else
              aR(i,N(ng))=2.0_r8*Bio(i,N(ng),ibio)-aL(i,N(ng))
# endif
              dR(i,N(ng))=(2.0_r8*aR(i,N(ng))+aL(i,N(ng))-              &
     &                     3.0_r8*Bio(i,N(ng),ibio))**2
              dL(i,N(ng))=(3.0_r8*Bio(i,N(ng),ibio)-                    &
     &                     2.0_r8*aL(i,N(ng))-aR(i,N(ng)))**2
!
              aR(i,1)=aL(i,2)
# ifdef NEUMANN
              aL(i,1)=1.5_r8*Bio(i,1,ibio)-0.5_r8*aR(i,1)
# else
              aL(i,1)=2.0_r8*Bio(i,1,ibio)-aR(i,1)
# endif
              dR(i,1)=(2.0_r8*aR(i,1)+aL(i,1)-                          &
     &                 3.0_r8*Bio(i,1,ibio))**2
              dL(i,1)=(3.0_r8*Bio(i,1,ibio)-                            &
     &                 2.0_r8*aL(i,1)-aR(i,1))**2
            END DO
#endif
!
! Reconcile interfacial values aR and aL using Weighted Essentially
! Non-Oscillatory (WENO) procedure.
!
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                deltaL=MAX(dL(i,k  ),eps)
                deltaR=MAX(dR(i,k+1),eps)
                r(i,k)=(deltaR*aR(i,k)+deltaL*aL(i,k+1))/               &
     &                 (deltaR+deltaL)
              END DO
            END DO
            DO i=Istr,Iend
#ifdef NEUMANN
              r(i,N(ng))=1.5_r8*Bio(i,N(ng),ibio)-0.5_r8*r(i,N(ng)-1)
              r(i,0    )=1.5_r8*Bio(i,1    ,ibio)-0.5_r8*r(i,1      )
#else
              r(i,N(ng))=2.0_r8*Bio(i,N(ng),ibio)-r(i,N(ng)-1)
              r(i,0    )=2.0_r8*Bio(i,1    ,ibio)-r(i,1      )
#endif
            END DO
!
! Remapping step: This operation consists essentially of three stages:
!---------------- (1) within each grid box compute averaged slope
! (stored as dR) and curvature (stored as dL); then (2) compute
! interfacial fluxes FC; and (3) apply these fluxes to complete
! remapping step.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
#ifdef LIMIT_INTERIOR
                deltaR=r(i,k)-Bio(i,k,ibio)      ! Constrain parabolic
                deltaL=Bio(i,k,ibio)-r(i,k-1)    ! segment monotonicity
                cffR=2.0_r8*deltaR               ! like in PPM
                cffL=2.0_r8*deltaL
                IF (deltaR*deltaL.lt.0.0_r8) THEN
                  deltaR=0.0_r8
                  deltaL=0.0_r8
                ELSE IF (ABS(deltaR).gt.ABS(cffL)) THEN
                  deltaR=cffL
                ELSE IF (ABS(deltaL).gt.ABS(cffR)) THEN
                  deltaL=cffR
                END IF
                aR(i,k)=Bio(i,k,ibio)+deltaR
                aL(i,k)=Bio(i,k,ibio)-deltaL
#else
                aR(i,k)=r(i,k  )
                aL(i,k)=r(i,k-1)
#endif
                dL(i,k)=0.5_r8*(aR(i,k)-aL(i,k))
                dR(i,k)=0.5_r8*(aR(i,k)+aL(i,k))-Bio(i,k,ibio)
              END DO
            END DO
!
!  Compute interfacial fluxes. The convention is that Wbio is positive
!  for upward motion (swimming, floating) and negative for downward motion
!  (sinking).
!
            wdt=-Wbio(iswim)*dtdays
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                IF (wdt.gt.0.0_r8) THEN           ! downward vertical
                  alpha=Hz(i,j,k+1)               ! motion (sinking)
                  cff =aL(i,k+1)
                  cffL=dL(i,k+1)
                  cffR=dR(i,k+1)
                  dz=wdt
                ELSE                              ! upward vertical
                  alpha=-Hz(i,j,k)                ! motion (swimming,
                  cff =aR(i,k)                    ! migration)
                  cffL=-dL(i,k)
                  cffR=dR(i,k)
                  dz=wdt
!!                IF (ABS(z_w(i,j,k)).lt.21.0_r8) THEN
!!                  dz=wdt*(1.0_r8-TANH((21.0_r8+z_w(i,j,k))*0.1_r8))
!!                ELSE
!!                  dz=wdt
!!                END IF
                END IF
                alpha=dz/alpha                    ! Courant number
                FC(i,k)=dz*(cff+alpha*(cffL-cffR*(3.0_r8-2.0_r8*alpha)))
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0    )=0.0_r8
              FC(i,N(ng))=0.0_r8
            END DO
!
!  Add semi-Lagrangian vertical flux.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                Bio(i,k,ibio)=Bio(i,k,ibio)+cff
              END DO
            END DO

          END DO SWIM_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentrations affected by BGC processes. This also takes
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

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
