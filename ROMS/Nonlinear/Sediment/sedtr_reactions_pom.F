#include "cppdefs.h"

      MODULE sedtr_reactions_pom_mod

#if defined NONLINEAR && defined SEDIMENT  && defined SEDTR_REACTIONS \
                     && defined SEDBIO_COUP
!
!svn $Id: sedtr_reactions_pom_bed.F 830 2017-01-24 21:21:11Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes reaction term for particles and particle-     !
!  associated material.                                                !
!                                                                      !
!  Module built based on version of Soetaert model as coded by K.      !
!  Fennel and Robin West (Dalhousie), and was incorporated into ROMS by!
!                                                                      !
!  Moriarty, J.M., C.K. Harris, K. Fennel, M.A.M. Friedrichs, K. Xu,   !
!    and C. Rabouille, 2017: The roles of resuspension, diffusion and  !
!    biogeochemical processes on oxygen dynamics offshore of the Rhone !
!    River, France: a numerical modeling study, Biogeosciences, 14,    !
!    1919-1946.                                                        !
!                                                                      !
!  Currently, reactions for particulate organic matter and dissolved   !
!  nutrients are available.  At present we include multiple organic    !
!  carbon fractions, oxygen, nitrate,   ammonium, and oxygen demand    !
!  units (ODUs).  We account for oxic and anerobic mineralization, and !
!  denitrification, as well as nitrification and oxidation of ODUs.    !
!  Annamox and phosphorus code could be added later.  Changes in       !
!  porosity are neglected. Inert organic matter should be modeled using!
!  inert sediment classes. Since publication, the model has been       !
!  adapted to allow particulate organic matter to affect seabed mass   !
!  and seabed layer thicknesses, but neglects variations in N:C ratios !
!  and biodiffusion coefficients.                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sedtr_reactions_pom

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sedtr_reactions_pom (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_biology
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
      CALL wclock_on (ng, iNLM, 16, __LINE__, __FILE__)
# endif
      CALL sedtr_reactions_pom_tile (ng, tile,                          &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   OCEAN(ng) % t,                                 &
     &                   SEDBED(ng) % bed_mass,                         &
     &                   SEDBED(ng) % bed_frac,                         &
     &                   SEDBED(ng) % bed)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE sedtr_reactions_pom
!
!***********************************************************************
      SUBROUTINE sedtr_reactions_pom_tile (ng, tile,                    &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         t,                                       &
     &                         bed_mass, bed_frac, bed)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
      USE mod_biology
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
# endif
!
!  Local variable declarations.
!
      integer :: i, ised, j, k

      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: zeros = 0.0_r8
      real(r8), parameter :: ones = 1.0_r8
      real(r8) :: cff

      real(r8) :: o2,no3, nh4, odu
      real(r8) :: oc(NST)
      real(r8) :: docdt(NST)
      real(r8) :: docdt_tot, doc1dt, doc2dt, do2dt, dodudt
      real(r8) :: dnh4dt, dno3dt
      real(r8) :: sConv, pwConv, poroConv
      real(r8) :: rxnOC1, rxnOC2
      real(r8) :: rxnOC1_o2, rxnOC2_o2, rxnOC_o2
      real(r8) :: rxnNH4_o2, rxnODU_o2, rxnOC_odu
      real(r8) :: rxnOC_no3, rxnNH4_amx
      real(r8) :: o2lim, totlim, odulim, dnflim
      real(r8) :: oxodu,nit
      real(r8) :: rem

!     Coefficients for remineralization and oxidation reactions
      real(r8) :: r_om=0.15_r8             !ratio mol N: mol C for POM
      real(r8) :: kg2mmolC = 1000.0_r8/12.0107_r8 !ratio of
                                           !sediment mass to mmol Carbon
                                           ![mmol/kg]
      real(r8) :: k_o2 = 1.0_r8            !Half-sat for O2 lim in 
                                           !aerobic remin
                                           ![umol o2/L == mmol/m3]
      real(r8) :: kin_odu = 1.0_r8         !Half-sat for O2 inhib. of 
                                           !anaerobic remin
                                           ![umol o2/L == mmol/m3]
      real(r8) :: oxodu_base = 20.0_r8     !oxidation rate of ODUs [1/d]
      real(r8) :: k_odu = 1.0_r8           !Half-sat for O2 lim in 
                                           !ODU oxidation
                                           ![umol o2/L == mmol/m3]
      real(r8) :: nit_base = 100.0_r8      !Nitrification Rate [1/d]
      real(r8) :: k_nit = 10.0_r8          !Half-sat for O2 lim of 
                                           !nitrification
                                           ![umol o2/L == mmol/m3]
      real(r8) :: k_dnf = 20.0_r8          !Half-sat for NO3 lim in 
                                           !denitrification
                                           ![umol o2/L == mmol/m3]
      real(r8) :: kin_dnf = 1.0_r8         !Half-sat for O2 inhib of 
                                           !denitrification
                                           ![umol o2/L == mmol/m3]
      real(r8) :: kin_anox_no3 = 10.0_r8   !Half-sat for NO3 inhib of 
                                           !anerobic remin
                                           ![umol no3/L == mmol/m3]
      real(r8) :: PB = 0.995_r8            !Fraction of newly created 
                                           !ODUS that are
                                           !assumed to be inert
      real(r8) :: NH4ads = 1.3_r8          !NH4 adsorption coefficient 
# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Particulate organic matter remineralization & 
!  Update to bed mass.
!-----------------------------------------------------------------------
!
      J_LOOP : DO j=Jstr,Jend
        I_LOOP : DO i=Istr,Iend
          K_LOOP : DO k=1,Nbed

!
!-----------------------------------------------------------------------
! Prepare for Soetaert routines
!-----------------------------------------------------------------------
!
! Conversion factors to account for porosity and bed thickness
! units are meters for sConv and pwConv
             sConv  = (1.0_r8-bed(i,j,k,iporo))*bed(i,j,k,ithck)
             pwConv = bed(i,j,k,iporo)*bed(i,j,k,ithck)
             poroConv = (1.0_r8-bed(i,j,k,iporo))/bed(i,j,k,iporo)

! Calculate rate constants for nitrification, ODU oxidation (1/s)
             nit=nit_base/24.0_r8/3600.0_r8 
             oxodu=oxodu_base/24.0_r8/3600.0_r8 

! Convert tracer units from ROMS to Soetaert units
! From mmol/m2 (bulk) to  mmol / m3 solids or porewater
! Also convert POM from kg to mmol C, and set oc = 0 for inorganics
!
!IF ((i.eq.3).and.(j.eq.3).and.(k.eq.1)) THEN
!write(6,*),"w.c.conc 1",t(i,j,1,nnew,iwc)/Hz(i,j,1)
!write(6,*),"sd conc 10",bed(i,j,1,iboxy)/                               &
!     &                  bed(i,j,1,ithck)/bed(i,j,1,iporo)
!ENDIF
             o2 =bed(i,j,k,iboxy)/pwConv 
             no3=bed(i,j,k,ibno3)/pwConv
             nh4=bed(i,j,k,ibnh4)/pwConv
             odu=bed(i,j,k,ibodu)/pwConv
!IF ((i.eq.3).and.(j.eq.3).and.(k.eq.1)) THEN
!write(6,*) "bed_oxy: ",o2
!write(6,*) "bed_no3: ",no3
!write(6,*) "bed_nh4: ",nh4
!write(6,*) "bed_odu: ",odu
!END IF
             DO ised=1,NST
               IF (sed_rxn(ised,ng).ne.0.0_r8) THEN
                 oc(ised)=bed_mass(i,j,k,nnew,ised)*kg2mmolC/sConv
               ELSE
                 oc(ised)=0.0_r8
               END IF
             END DO

! Initialize derivatives
             dnh4dt = 0.0_r8
             dno3dt = 0.0_r8
             do2dt  = 0.0_r8
             dodudt = 0.0_r8
             docdt_tot = 0.0_r8

! Calculate limitation terms that are used to partition OM
! remineralization between aerobic, nitrate and anerobic 
! remineralization
             o2lim = o2/(o2+k_o2)
             odulim = kin_odu/(o2+kin_odu)                              &
     &            *kin_anox_NO3/(no3+kin_anox_NO3)
             dnflim = no3/(no3+k_dnf)                                   &
     &            *kin_dnf/(kin_dnf+o2)
             totlim = o2lim+odulim+dnflim+eps

! Calculate rate constant for POM remineralization (1/s),
! change in POM due to remineralization,
! and change to bed mass due to decaying POM
             DO ised=1,NST
              rem=(sed_rxn(ised,ng)                                     &
# ifdef SEDBIO_POM_TVAR
     &            + 1.0_r8**(t(i,j,1,nstp,itemp)-20.0_r8)               &
# endif
     &            )/24.0_r8/3600.0_r8
              docdt(ised) = 0.0_r8 + oc(ised)*rem
              docdt_tot = docdt_tot + docdt(ised)
            END DO

! Calculate O2 consumed by aerobic respiration (O2:C = 1:1)
! in mmol O2/m3/s porewater
            rxnOC_o2  = docdt_tot*o2lim/totlim *poroConv

! Calculate NO3 consumed by denitrification (NO3:C = 0.8:1)
! (mmol NO3/m3 porewater /sec
            rxnOC_no3 = 0.8_r8*docdt_tot*dnflim/totlim*poroConv

! Calculate ODUs produced by anaerobic metabolism (ODU:C = 1:1).
! (mmol ODU/m3 porewater/sec)
! Note that a fraction (=PB) of the ODUs are considered to be lost
! because they are assumed to be non-reactive solids.
            rxnOC_odu = docdt_tot*odulim/totlim                         &
     &            *poroConv*(1.0_r8-PB)

! Calculate O2 consumed by nitrification (O2:NH3 = 2:1) (rxnNH4_o2)
! and by ODU oxidation (O2:ODU = 1:1) (rxnODU_o2)
! in mmol O2/m3 porewater /sec
!
! Note that calculations assume that remineralizaiton occurs first
! so rxnOC_* are subtracted/added to o2, no3, and odu
            rxnNH4_o2 = 2.0_r8*nit*(nh4+docdt_tot*r_om*poroConv)        &
     &            *(o2-rxnOC_o2)/((o2-rxnOC_o2)+k_nit)
!           rxnNH4_o2 = 2.0_r8*nit*(no3-rxnOC_no3)*(o2-rxnOC_o2)        &
!    &            /((o2-rxnOC_o2)+k_nit)
!           rxnNH4_o2 = 2.0_r8*nit*nh4*o2/(o2+k_nit)
            rxnODU_o2 = oxodu*(o2-rxnOC_o2)*(odu+rxnOC_odu)             &
     &            /((o2-rxnOC_o2)+k_odu)
!           rxnODU_o2 = oxodu*o2*odu/(o2+k_odu)

! Calculate NH4 consumed by anammox: 5NH3 + 3HNO3 -> 4N2 + 9H2O
! simplify: Neglect anammox
            rxnNH4_amx = 0.0_r8

!
!-----------------------------------------------------------------------
! Update derivatives of dissolved tracers
!-----------------------------------------------------------------------
! All units are mmol tracer /m3 porewater/s
!
             do2dt = do2dt - rxnOC_o2 - rxnNH4_o2 - rxnODU_o2
             dodudt = dodudt + rxnOC_odu - rxnODU_o2
             dno3dt = dno3dt - rxnOC_no3 + rxnNH4_o2/2.0_r8             &
     &            - rxnNH4_amx*3.0_r8/5.0_r8
             dnh4dt = dnh4dt - rxnNH4_o2/2.0_r8 - rxnNH4_amx +          &
     &            docdt_tot*r_om*poroConv
! Ammonium adsorption
             dnh4dt=dnh4dt/(1.0_r8+NH4ads)

! Use updated derivatives to update tracer values.
! Convert from mmol/m3 porewater or solids to mmol / m2 (bulk)
             o2  = (o2  + dt(ng)*do2dt) *pwConv
             no3 = (no3 + dt(ng)*dno3dt)*pwConv
             nh4 = (nh4 + dt(ng)*dnh4dt)*pwConv
             odu = (odu + dt(ng)*dodudt)*pwConv
             DO ised=1,NST
               oc(ised) = (oc(ised) - dt(ng)*docdt(ised))*sConv         &
     &               /kg2mmolC
             END DO

! Make sure tracers are positive
             o2=MAX(o2,0.0_r8)
             no3=MAX(no3,0.0_r8)
             nh4=MAX(nh4,0.0_r8)
             odu=MAX(odu,0.0_r8)
             DO ised=1,NST
                oc(ised)=MAX(oc(ised),0.0_r8)
             END DO
!
! Map bed tracer values back into ROMS variable
!
             bed(i,j,k,iboxy)=o2
             bed(i,j,k,ibno3)=no3
             bed(i,j,k,ibnh4)=nh4
             bed(i,j,k,ibodu)=odu
             DO ised=1,NST
               IF (sed_rxn(ised,ng).ne.0.0_r8) THEN
                 bed_mass(i,j,k,nnew,ised)=oc(ised)
               ELSE
                 bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nnew,ised)
               END IF
             END DO
!IF ((i.eq.3).and.(j.eq.3).and.(k.eq.1)) THEN
!write(6,*),"w.c.conc 1",t(i,j,1,nnew,iwc)/Hz(i,j,1)
!write(6,*),"sd conc 15",bed(i,j,1,iboxy)/                               &
!     &                  bed(i,j,1,ithck)/bed(i,j,1,iporo)
!ENDIF

          END DO K_LOOP
        END DO I_LOOP
      END DO J_LOOP


!
!-----------------------------------------------------------------------
!  Update bed properties following reaction terms.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,Nbed
            bed(i,j,k,ithck)=zeros
            cff=zeros
            DO ised=1,NST
              cff=cff+bed_mass(i,j,k,nnew,ised)
            END DO
            cff=MAX(eps,cff)
            DO ised=1,NST
              bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff
              bed(i,j,k,ithck)=bed(i,j,k,ithck)+                        &
     &            bed_mass(i,j,k,nnew,ised)/(Srho(ised,ng)*             &
     &            (ones-bed(i,j,k,iporo)))  
            END DO
!IF ((i.eq.3).and.(j.eq.3).and.(k.eq.1)) THEN
!write(6,*),"w.c.conc 1",t(i,j,1,nnew,iwc)/Hz(i,j,1)
!write(6,*),"sd conc 19",bed(i,j,1,iboxy)/                               &
!     &                  bed(i,j,1,ithck)/bed(i,j,1,iporo)
!ENDIF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
# endif

      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed)
# endif
      RETURN
      END SUBROUTINE sedtr_reactions_pom_tile
#endif
      END MODULE sedtr_reactions_pom_mod
