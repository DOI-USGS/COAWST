#include "cppdefs.h"

      MODULE sedtr_reactions_sed_decay_mod

#if defined NONLINEAR && defined SEDIMENT  && defined SEDTR_REACTIONS \
                     && defined SEDTR_SED_DECAY
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
      PUBLIC  :: sedtr_reactions_sed_decay

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sedtr_reactions_sed_decay (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
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
      CALL sedtr_reactions_sed_decay_tile (ng, tile,                    &
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
      END SUBROUTINE sedtr_reactions_sed_decay
!
!***********************************************************************
      SUBROUTINE sedtr_reactions_sed_decay_tile (ng, tile,              &
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

      real(r8) :: oc(NST)
      real(r8) :: docdt(NST)
      real(r8) :: rem

!     Coefficients for remineralization and oxidation reactions
      real(r8) :: r_om=0.15_r8             !ratio mol N: mol C for POM
      real(r8) :: kg2mmolC = 1000.0_r8/12.0107_r8 !ratio of
                                           !sediment mass to mmol Carbon
                                           ![mmol/kg]
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

! Convert tracer units from ROMS to Soetaert units
! From mmol/m2 (bulk) to  mmol / m3 solids or porewater
! Also convert POM from kg to mmol C, and set oc = 0 for inorganics
!
            DO ised=1,NST

              IF (sed_rxn(ised,ng).ne.0.0_r8) THEN
! Calculate rate constant for POM remineralization (1/s),
! change in POM due to remineralization,
! and change to bed mass due to decaying POM 
                rem=(sed_rxn(ised,ng)                                   &
# ifdef SEDBIO_POM_TVAR
     &            + 1.0_r8**(t(i,j,1,nstp,itemp)-20.0_r8)               &
# endif
     &            )/24.0_r8/3600.0_r8
                docdt(ised) = 0.0_r8 + bed_mass(i,j,k,nnew,ised)*rem
! Use updated derivatives to update tracer values.
! Convert from mmol/m3 porewater or solids to mmol / m2 (bulk)
                bed_mass(i,j,k,nnew,ised) = MAX(0.0_r8,                 &
     &             bed_mass(i,j,k,nnew,ised) - dt(ng)*docdt(ised))
                t(i,j,k,nnew,ised) = MAX(0.0_r8,                 &
     &             t(i,j,k,nnew,ised)*(1.0_r8 - dt(ng)*rem))



!              ELSE
!               oc(ised)=0.0_r8
!               bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nnew,ised)
              END IF

            END DO

!
!-----------------------------------------------------------------------
!  Update bed properties following reaction terms.
!-----------------------------------------------------------------------
!
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

          END DO K_LOOP
        END DO I_LOOP
      END DO J_LOOP
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
      END SUBROUTINE sedtr_reactions_sed_decay_tile
#endif
      END MODULE sedtr_reactions_sed_decay_mod
