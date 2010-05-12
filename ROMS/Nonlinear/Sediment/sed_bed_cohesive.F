#include "cppdefs.h"

      MODULE bed_cohesive_mod
#if defined NONLINEAR && defined SEDIMENT && defined COHESIVE_BED
!
!svn $Id: sediment.F 1577 2008-07-22 19:16:57Z bferre $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  the following:                                                      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller.                             !
!  * Bedload based on Soulsby combined waves + currents                !
!    (p166 Soulsby 1997)                                               !
!  * Bedload slope term options: Nemeth et al, 2006, Coastal           !
!    Engineering, v 53, p 265-275; Lesser et al, 2004, Coastal         !
!    Engineering, v 51, p 883-915.                                     !
!                                                                      !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: bed_cohesive

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bed_cohesive (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
# ifdef BBL_MODEL
      USE mod_bbl
# endif
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
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL bed_cohesive_tile (ng, tile,                                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng), nnew(ng),                           &
# ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
# endif
# ifdef BBL_MODEL
     &                    BBL(ng) % bustrc,                             &
     &                    BBL(ng) % bvstrc,                             &
     &                    BBL(ng) % bustrw,                             &
     &                    BBL(ng) % bvstrw,                             &
     &                    BBL(ng) % bustrcwmax,                         &
     &                    BBL(ng) % bvstrcwmax,                         &
# endif
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
# if defined SED_MORPH
     &                    GRID(ng) % bed_thick,                         &
# endif
     &                    OCEAN(ng) % bed,                              &
     &                    OCEAN(ng) % bed_frac,                         &
     &                    OCEAN(ng) % bed_mass,                         &
     &                    OCEAN(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE bed_cohesive
!
!***********************************************************************
      SUBROUTINE bed_cohesive_tile (ng, tile,                           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp, nnew,                             &
# ifdef WET_DRY
     &                          rmask_wet,                              &
# endif
# ifdef BBL_MODEL
     &                          bustrc, bvstrc,                         &
     &                          bustrw, bvstrw,                         &
     &                          bustrcwmax, bvstrcwmax,                 &
# endif
     &                          bustr, bvstr,                           &
# if defined SED_MORPH
     &                          bed_thick,                              &
# endif
     &                          bed, bed_frac, bed_mass,                &
     &                          bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef BEDLOAD
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
      real(r8), intent(in) :: bustrw(LBi:,LBj:)
      real(r8), intent(in) :: bvstrw(LBi:,LBj:)
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: Ksed, i, indx, ised, j, k, ks
      integer :: bnew

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      real(r8) :: cu, cff4
      real(r8) :: thck_avail, thck_to_add

      integer, dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: ksource

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,NST) :: dep_mass
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: tau_w
# if defined SED_MORPH
      integer :: kbed
# endif
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Qb
# if defined MIXED_BED
      real(r8) :: pcoh 
# endif
# if defined COHESIVE_BED || defined MIXED_BED
      real(r8) :: alpha, bzactv, frt, tcb_temp, tcb_top, tcb_bot
      real(r8), dimension(Nbed) :: tcr
# endif
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      real(r8), dimension(Nbed) :: bmz
      real(r8), dimension(Nbed) :: bz
# endif

# include "set_bounds.h"

# ifdef BEDLOAD
      bnew=nnew
# else
      bnew=nstp
# endif
!-----------------------------------------------------------------------
! Compute maximum bottom stress.
!-----------------------------------------------------------------------
!
# if defined BEDLOAD_MPM || defined SUSPLOAD
#  ifdef BBL_MODEL
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrcwmax(i,j)*bustrcwmax(i,j)+              &
     &                    bvstrcwmax(i,j)*bvstrcwmax(i,j))
#   ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#   endif
        END DO
      END DO
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE Istr-1,Iend+1
#   else
#    define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
#   ifdef NS_PERIODIC
#    define J_RANGE Jstr-1,Jend+1
#   else
#    define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
      DO i=I_RANGE
        DO j=J_RANGE
          tau_w(i,j)=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*             &
     &                           (bustr(i,j)+bustr(i+1,j))+             &
     &                           (bvstr(i,j)+bvstr(i,j+1))*             &
     &                           (bvstr(i,j)+bvstr(i,j+1)))
#   ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#   endif
        END DO
      END DO
#   undef I_RANGE
#   undef J_RANGE
#  endif
# endif
!
# ifdef SUSPLOAD
!
      J_LOOP : DO j=Jstr,Jend
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
!
!-----------------------------------------------------------------------
!  STEP 1 Sediment deposition and resuspension near the bottom.
!-----------------------------------------------------------------------
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precipitation flux FC(:,0), already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
          DO i=Istr,Iend
#  ifdef SED_MORPH
!
! Apply morphology factor.
!
            ero_flux(i,j,ised)=ero_flux(i,j,ised)*morph_fac(ised,ng)
            settling_flux(i,j,indx)=settling_flux(i,j,indx)*morph_fac(ised,ng)
#  endif /* defined SED_MORPH */
!
!  Record mass deposited (or eroded) for evaluation later
!
            dep_mass(i,ised)=-(ero_flux(i,ised)-settling_flux(i,j,ised))
!
!  Add or subtract material to top layer and update bed mass arrays
!
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,bnew,ised)+    &
     &                             dep_mass(i,ised),0.0_r8)
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
          END DO
        END DO SED_LOOP
!
!-----------------------------------------------------------------------
!  STEP 2. At this point, all deposition or erosion is complete, and
!  has been added/subtracted to top layer. Thickness has NOT been corrected.
!-----------------------------------------------------------------------
!
        DO i=Istr,Iend
!
!  Determine if *net* deposition ocurred
!
          cff=0.0_r8
          cff2=0.0_r8
          DO ised=1,NST
             cff=cff+dep_mass(i,ised)
             cff2=cff2+ero_flux(i,ised)
          END DO
          bottom(i,j,idnet)=cff
          bottom(i,j,idflux)=cff2
          IF ( bottom(i,j,idnet).GT.0.0_r8) THEN

#if defined COHESIVE_BED || defined MIXED_BED
!
!  Calculate tau_crit of deposited bed
!
!  Calculate new mass in top layer
             bmz(1) = 0.0_r8
             bmz(2) = 0.0_r8
             DO ised=1,NST
                bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
                bmz(2) = bmz(2)+bed_mass(i,j,2,nnew,ised)
             END DO
#  undef JP_TCR
#  undef BF_TCR
#  if defined JP_TCR
!  Set top layer critical shear stress to
!  weighted average of tcr_min and linear profile in first layer 
             cff1 = bottom(i,j,idnet)/MAX(bmz(1),eps) 
             cff2 = bed(i,j,1,ibtcr)+bed(i,j,2,ibtcr)-2*tcr_min(ng)
             bed(i,j,1,ibtcr) = bed(i,j,1,ibtcr) - (cff1*cff2)
#  elif defined BF_TCR
             cff = 0.25
             bed(i,j,1,ibtcr) = (1-cff)*rho0*tau_w(i,j)+cff*tcr_min(ng)
#  else
             IF (Nbed.GT.1) THEN
!  Average of cff1 and cff2, where
!    cff1 = linear extension of previous tcr slope to new surface
!    cff2 = minimum deposition
                cff1 = bed(i,j,1,ibtcr) -                               &
     &             bottom(i,j,idnet) *                                  &
     &             (bed(i,j,2,ibtcr)-bed(i,j,1,ibtcr)) /                &
     &             (bmz(1)-bottom(i,j,idnet))
                cff2 = MAX( rho0*tau_w(i,j) , tcr_min(ng) )
                bed(i,j,1,ibtcr) = MIN(bed(i,j,1,ibtcr),                &
     &                                 MAX( 0.5*(cff1+cff2), cff2 ))
             ELSE
!  TODO: Not sure what mud tau_crit should be for single-layer bed.
!  Try weighted average of dep and old value
                cff1 = (bed(i,j,1,ibtcr)*(bmz(1)-bottom(i,j,idnet))     &
     &                  + cff2*bottom(i,j,idnet)) / bmz(1)
                bed(i,j,1,ibtcr) = MAX( cff1, cff2 )
             END IF
#  endif
#endif
!
!  Deposition. Determine if we need to create a new bed layer 
! TODO: fix kludge...this will break if NBED < 3
            IF ((time(ng).gt.(bed(i,j,1,iaged)+1.1_r8*dt(ng))).and.     &
     &         (bed(i,j,1,ithck).gt.newlayer_thick(ng)))                &
     &        THEN

!
!  TODO - Age of bed should be weighted by sediment mass...but this will
!  mess up the test above
              bed(i,j,1,iaged)=time(ng)
              IF (Nbed.gt.1) THEN
                 IF(bed(i,j,2,ithck).lt.minlayer_thick(ng)) THEN
!
! Layer 2 is smaller than minimum size
! Instead of pushing down all layers, just combine top layers
!

                    cff=0.0_r8
                    cff1=0.0_r8
                    cff2=0.0_r8
                    DO ised=1,NST
                       cff =cff +dep_mass(i,ised)
                       cff1=cff1+bed_mass(i,j,1,nnew,ised)
                       cff2=cff2+bed_mass(i,j,2,nnew,ised)
                    END DO

! Combine poro & age
! TODO: Need to update porosity and age
#if defined COHESIVE_BED || defined MIXED_BED
!
!   Assign new tau_crit at 2nd layer
!
!                    bed(i,j,2,ibtcr)= bed(i,j,2,ibtcr) -                &
!     &       (cff1-cff)/(cff1+cff2)*(bed(i,j,3,ibtcr)-bed(i,j,1,ibtcr))
                  bed(i,j,2,ibtcr) = 0.5_r8*( bed(i,j,1,ibtcr)+         &
     &                                        bed(i,j,2,ibtcr) )
#endif
!
!  Update bed mass
!
          DO ised=1,NST
             bed_mass(i,j,2,nnew,ised)=MAX(bed_mass(i,j,2,nnew,ised)+   &
     &            bed_mass(i,j,1,nnew,ised)-dep_mass(i,ised),0.0_r8)
             bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
          END DO 

       ELSE
!
!  Combine bottom layers.
              bed(i,j,Nbed,iporo)=0.5_r8*(bed(i,j,Nbed-1,iporo)+        &
     &                                    bed(i,j,Nbed,iporo))
              bed(i,j,Nbed,iaged)=0.5_r8*(bed(i,j,Nbed-1,iaged)+        &
     &                                    bed(i,j,Nbed,iaged))
              DO ised=1,NST
                bed_mass(i,j,Nbed,nnew,ised)=                           &
     &                             bed_mass(i,j,Nbed-1,nnew,ised)+      &
     &                             bed_mass(i,j,Nbed  ,nnew,ised)
              END DO
#if defined COHESIVE_BED || defined MIXED_BED
                    cff1=0
                    cff2=0
                    DO ised=1,NST
                       cff1=cff1+bed_mass(i,j,Nbed-1,nnew,ised)
                       cff2=cff2+bed_mass(i,j,Nbed,nnew,ised)
                    END DO
!
!   Assign tcrit at top of new bottom bed tcrit for Nbed-1
              bed(i,j,Nbed,ibtcr)= bed(i,j,Nbed-1,ibtcr)
!              bed(i,j,Nbed,ibtcr)=bed(i,j,Nbed-1,ibtcr) +                &
!     &            cff2*(bed(i,j,Nbed,ibtcr)-bed(i,j,Nbed-1,ibtcr))/cff1
#endif
!
!  Push layers down.
!
              DO k=Nbed-1,2,-1
                bed(i,j,k,iporo)=bed(i,j,k-1,iporo)
                bed(i,j,k,iaged)=bed(i,j,k-1,iaged)
                DO ised =1,NST
                  bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k-1,nnew,ised)
                END DO
#if defined COHESIVE_BED || defined MIXED_BED
                bed(i,j,k,ibtcr)=bed(i,j,k-1,ibtcr)
#endif
              END DO    !k=Nbed-1,2,-1
!
!  Set new top layer parameters.
!
#if defined COHESIVE_BED || defined MIXED_BED
!              bed(i,j,1,ibtcr)=tcb_temp
!#  if defined JP_TCR
!
!  Tau_crit at top already determined
!  Set tau_crit at 2nd layer
!
          cff=0.0_r8
          cff1=0.0_r8
          DO ised=1,NST
             cff =cff +dep_mass(i,ised)
             cff1=cff1+bed_mass(i,j,1,nnew,ised)
          END DO

          cff2 = (bed(i,j,2,ibtcr)-bed(i,j,1,ibtcr))/cff1  !BF corrected
          bed(i,j,2,ibtcr) = bed(i,j,1,ibtcr)+cff*cff2
!#  endif
#endif
              DO ised=1,NST
                bed_mass(i,j,2,nnew,ised)=MAX(bed_mass(i,j,2,nnew,ised)-&
     &                                    dep_mass(i,ised),0.0_r8)
                bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
              END DO
                 END IF    !bed(i,j,2,ithck).lt.minlayer_thick
            END IF    !NBED=1
          ELSE          ! No new bed layer
!
! Deposition has occured, but no new bed layer was created 
!            
#if defined COHESIVE_BED || defined MIXED_BED
! Set top layer critical shear stress to tcb_temp
!             bed(i,j,1,ibtcr)=tcb_temp
!             write(*,*) 'after dep, nnb: ',bed(i,j,1,ibtcr)
#endif
          END IF    ! New bed layer
        ELSEIF (bottom(i,j,idnet).LT.0.0_r8) THEN
!
!   Net Erosion Occurred
!
#if defined COHESIVE_BED || defined MIXED_BED
              bmz(1) = 0.0_r8
              DO ised=1,NST
                 bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
              END DO
!  recalc tc for top of bed, based on linear
!  interpolation and mass removed / orig. mass in top layer
!              write(*,*) 'before ero: ',bed(i,j,1,ibtcr),' cff: ',cff
              bed(i,j,1,ibtcr)=bed(i,j,1,ibtcr)+                        &
     &             MIN(1.0_r8,-bottom(i,j,idnet)/MAX(eps,bmz(1)))*      &
     &             MAX(0.0_r8,(bed(i,j,2,ibtcr)-bed(i,j,1,ibtcr)))     
!              write(*,*) 'after ero: ',bed(i,j,1,ibtcr)
#endif
          END IF    !cff.LT.0.0_r8

!
! Recalculate thickness and fractions for all layers.
!
          DO k=1,Nbed
            cff3=0.0_r8
            DO ised=1,NST
              cff3=cff3+bed_mass(i,j,k,nnew,ised)
            END DO
            IF (cff3.eq.0.0_r8) THEN 
               cff3=eps 
            END IF 
            bed(i,j,k,ithck)=0.0_r8
            DO ised=1,NST
              bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff3
              bed(i,j,k,ithck)=MAX(bed(i,j,k,ithck)+                    &
     &                         bed_mass(i,j,k,nnew,ised)/(              &
     &                           Srho(ised,ng)*                         &
     &                          (1.0_r8-bed(i,j,k,iporo))),0.0_r8)
            END DO
          END DO
        END DO
      END DO J_LOOP
!
!  End of Sediment only section.
!
# endif /* defined Sediment only section */

      J_LOOP2 : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  STEP 4. Ensure top bed layer thickness is greater or equal than 
!  active layer thickness. If need to add sed to top layer, then entrain
!  from lower levels. Create new layers at bottom to maintain Nbed.
!-----------------------------------------------------------------------
!
        I_LOOP_CB: DO i=Istr,Iend
#undef JP_ACTIVE
#undef BF_ACTIVE
#if defined JP_ACTIVE
          bottom(i,j,iactv)= 6.0_r8*bottom(i,j,isd50)
# elif defined BF_ACTIVE_MUD
          bottom(i,j,iactv)= MAX( 0.0_r8,                               &
     &                          0.007_r8*                               &
     &                         (rho0*tau_w(i,j)-bottom(i,j,itauc)))+    &
     &                          6.0_r8*bottom(i,j,isd50)
# elif defined COHESIVE_BED
! Calculate active layer for muddy beds
          bottom(i,j,iactv)=MAX( 6.0_r8*bottom(i,j,isd50),              &
     &                          0.007_r8*                               &
     &                         (rho0*tau_w(i,j)-bottom(i,j,itauc))+     &
     &                          6.0_r8*bottom(i,j,isd50) )
# elif defined MIXED_BED
! Calculate active layer for mixed beds
          cff1 = MAX( 0.0_r8,                                           &
     &                0.007_r8*                                         &
     &               (rho0*tau_w(i,j)-bottom(i,j,itauc))) +             &
     &                6.0_r8*bottom(i,j,isd50)
          cff2 = 6.0_r8*bottom(i,j,isd50)
          bottom(i,j,iactv)= MAX( 6.0_r8*bottom(i,j,isd50),             &
     &         bottom(i,j,idprp)*cff1 + (1-bottom(i,j,idprp))*cff2 )
# else
! TODO: Correct non-cohesive active layer thickness for ripples
!          bottom(i,j,iactv)=MAX( 6.0_r8*bottom(i,j,isd50),              &
!     &                          0.007_r8*                               &
!     &                         (rho0*tau_w(i,j)-bottom(i,j,itauc))+     &
!     &                          6.0_r8*bottom(i,j,isd50) )
          bottom(i,j,iactv)= MAX( 6.0_r8*bottom(i,j,isd50),             &
     &                          0.007_r8*                               &
     &                         (rho0*tau_w(i,j)-bottom(i,j,itauc)))
# endif

#  if defined COHESIVE_BED || defined MIXED_BED
!  Find first layer with tc > tb
!  Remember, the critical stresses apply to the TOP of each layer
          Ksed = 1
          bzactv = 0.0_r8
          frt = 0.0_r8 ! CRS
          cff1 = rho0*tau_w(i,j)
          tcb_top = bed(i,j,1,ibtcr)
          tcb_bot = bed(i,j,2,ibtcr)
# if defined MIXED_BED
!         Calc. tau crit for bottom of next layer
!         Update cohesive property of seds in top layer
          cff3 = 0.0_r8
          cff4 = 1.0_r8
          DO ised=1,NCS
             cff3=cff3+bed_frac(i,j,1,ised)
             cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          DO ised=NCS+1,NST
             cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
          ENDDO   
          pcoh=min(max((cff3-transN(ng))/(transC(ng)-transN(ng))        &
     &      ,0.0_r8),1.0_r8)
          tcb_top = pcoh*bed(i,j,1,ibtcr)+(pcoh-1.0_r8)*cff4*rho0   !BF
# endif
          IF(cff1 .GT. tcb_top)THEN
!            Calculate tcb_temp for next layer
             tcb_bot = bed(i,j,Ksed+1,ibtcr)
# if defined MIXED_BED
!            Recalculate cohesive fraction and mean grain tau crit
!            Note that Ksed is used for grain props, and Ksed+1 for
!            bed tau crit at bottom of layer Ksed
             cff3 = 0.0_r8
             cff4 = 1.0_r8
             DO ised=1,NCS
                cff3=cff3+bed_frac(i,j,Ksed,ised)
                cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,Ksed,ised)
             END DO
             DO ised=NCS+1,NST
                cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,Ksed,ised)
             ENDDO   
!            Calculate cohesive behavior and blended tau crit
             pcoh=min(max((cff3-transN(ng))/(transC(ng)-transN(ng)),    &
     &            0.0_r8),1.0_r8)
             tcb_bot =pcoh*bed(i,j,Ksed+1,ibtcr)+(pcoh-1.0_r8)*cff4*rho0  !BF
# endif
             DO WHILE ( (Ksed.LE.(Nbed-1)) .AND.                        &
     &            (cff1 .GT. tcb_bot))
!               Add entire layer
                bzactv = bzactv + bed(i,j,Ksed,ithck)
                tcb_top = tcb_bot
                tcb_bot = bed(i,j,Ksed+1,ibtcr)
# if defined MIXED_BED
!               Recalculate cohesive fraction and mean grain tau crit
                cff3 = 0.0_r8
                cff4 = 1.0_r8
                DO ised=1,NCS
                   cff3=cff3+bed_frac(i,j,Ksed,ised)
                   cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,Ksed,ised)
                END DO
                DO ised=NCS+1,NST
                   cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,Ksed,ised)
                ENDDO
!               Calculate cohesive behavior and blended tau crit
                pcoh=min(max((cff3-transN(ng))/(transC(ng)-transN(ng)), &
     &               0.0_r8),1.0_r8)
                tcb_bot = pcoh*bed(i,j,Ksed+1,ibtcr)+(pcoh-1.0_r8)*cff4 &
     &                         *rho0                                      !BF
#endif
                Ksed = Ksed+1
             ENDDO
             frt = MAX(0.0_r8,(cff1-tcb_top) ) /                        &
     &            MAX( eps,                                             &
     &            ( tcb_bot-tcb_top ))
          ENDIF
          bzactv = bzactv+frt*bed(i,j,Ksed,ithck)
          bzactv = MAX( bzactv, 6.0_r8*bottom(i,j,isd50) ) !CRS
          bottom(i,j,iactv)=min(bottom(i,j,iactv),bzactv) !?CRS
!          bottom(i,j,iactv)=max(bottom(i,j,iactv),bzactv) !?CRS
#  endif /* defined COHESIVE_BED || defined MIXED_BED */

# ifdef SED_MORPH
!
! Apply morphology factor.
          bottom(i,j,iactv)=MAX(bottom(i,j,iactv)*morph_fac(1,ng),      &
     &                          bottom(i,j,iactv))
# endif
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                        ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                END IF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0_r8)/                  &
     &             MAX(bed(i,j,Ksed,ithck),eps)	        
              DO ised=1,NST
                cff1=0.0_r8
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,nnew,ised)
                END DO
                bed_mass(i,j,1   ,nnew,ised)=cff1-                      &
     &                                bed_mass(i,j,Ksed,nnew,ised)*cff2
                bed_mass(i,j,Ksed,nnew,ised)=                           &
     &                                bed_mass(i,j,Ksed,nnew,ised)*cff2
              END DO
!
!  Update thickness of fractional layer ksource_sed.
              bed(i,j,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_r8)

#if defined COHESIVE_BED || defined MIXED_BED
!
! Update tau_cr of fractional layer
	bed(i,j,Ksed,ibtcr) = bed(i,j,Ksed+1,ibtcr)-                    &
     &         cff2*(bed(i,j,Ksed+1,ibtcr)-bed(i,j,Ksed,ibtcr))
#endif

!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
              END DO
              IF (cff3.eq.0.0_r8) THEN
                 cff3=eps
              END IF
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
              END DO
!
!  Upate bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)

#  if defined COHESIVE_BED || defined MIXED_BED
                bed(i,j,k-ks,ibtcr)=bed(i,j,k,ibtcr)
#  endif

                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,nnew,ised)=bed_mass(i,j,k,nnew,ised)
                END DO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
# if defined COHESIVE_BED || defined MIXED_BED
#define BF_TCR
#if defined BF_TCR
          bmz(1) = 0.0_r8
          DO ised=1,NST
             bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
          ENDDO
          DO k=2,Nbed
             bmz(k) = bmz(k-1)
             DO ised=1,NST
                bmz(k)=bmz(k)+bed_mass(i,j,k,nnew,ised)
             ENDDO
          ENDDO
          tcr(1) = tcr_min(ng)
          DO k=2,Nbed
             tcr(k) = tcr_min(ng)
             IF (bmz(k-1).GT.eps) THEN
                tcr(k) = exp((log(bmz(k-1))-                            &
     &               bottom(i,j,idoff))/                                &
     &               bottom(i,j,idslp))
             ENDIF
             tcr(k) = MIN( MAX( tcr(k), tcr_min(ng)), tcr_max(ng) )
          ENDDO
#else
              cff1 = cff*(tcr_max(ng)-bed(i,j,Nbed-ks,ibtcr))
#endif
# endif
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)

#  if defined COHESIVE_BED || defined MIXED_BED
!                bed(i,j,k,ibtcr)=bed(i,j,Nbed-ks,ibtcr)
! Interpolate bottom layer to tau_crit_max
! TODO: should be the reference profile and not tau_crit_max
#if defined BF_TCR
                bed(i,j,k,ibtcr)=tcr(k);
#else
                bed(i,j,k,ibtcr)=bed(i,j,Nbed-ks,ibtcr)+                &
     &                  REAL(k-Nbed+ks,r8) * cff1
#endif
#  endif

                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,nnew,ised)=                            &
     &                             bed_mass(i,j,Nbed-ks,nnew,ised)*cff
                END DO
              END DO
            END IF  ! Nbed > 1
          END IF  ! increase top bed layer

#if defined MIXED_BED
!
!         Update cohesive property of seds in top layer
          cff1 = 0.0_r8
          DO ised=1,NCS
             cff1=cff1+bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,idprp)=min(max((cff1-transN(ng))/                  &
     &     (transC(ng)-transN(ng)),0.0_r8),1.0_r8)
#endif
        END DO I_LOOP_CB

#if defined COHESIVE_BED || defined MIXED_BED
         I_LOOP_CB2: DO i=Istr,Iend
!
!  Key to this algorithm: "mass available depth" (mad) [kg/m2]
!  present tau crit profile (tcp)
!  representative tau crit profile (tcr)          
!  Compute depth and mass depth of sediment column
!  Compute cumulative depth
!  Compute mass depth
          bmz(1) = 0.0_r8
          DO ised=1,NST
             bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
          ENDDO
          DO k=2,Nbed
             bmz(k) = bmz(k-1)
             DO ised=1,NST
                bmz(k)=bmz(k)+bed_mass(i,j,k,nnew,ised)
             ENDDO
          ENDDO

!  Calculate representative critical shear stress profile
!  Note that the values are for the TOP of the layer...we
!  assume the bottom of the bottom layer has tcr = tcr_max
          tcr(1) = tcr_min(ng)
          DO k=2,Nbed
             tcr(k) = tcr_min(ng)
             IF (bmz(k-1).GT.eps) THEN
                tcr(k) = exp((log(bmz(k-1))-                            &
     &               bottom(i,j,idoff))/                                &
     &               bottom(i,j,idslp))
             ENDIF
             tcr(k) = MIN( MAX( tcr(k), tcr_min(ng)), tcr_max(ng) )
          ENDDO
          DO k=Nbed,Nbed-ks,-1
             bed(i,j,k,ibtcr)=tcr(k);
          ENDDO
!  Relax tau crit profile bottom(i,j,k,ibtcr) 
!  towards representative profile tcr...100 x slower for "swelling"
!  than consolidation
!  TODO - make the factor an input parameter
          IF( bed(i,j,1,ibtcr).LE.tcr(1)) THEN
             alpha = MIN(dt(ng)/bottom(i,j,idtim),1.0_r8)
          ELSE
             alpha = MIN(dt(ng)/(100.0_r8*bottom(i,j,idtim)),1.0_r8)
          ENDIF
          DO k=1,Nbed
             bed(i,j,k,ibtcr)=bed(i,j,k,ibtcr)+                         &
     &            alpha*(tcr(k)-bed(i,j,k,ibtcr))
             bed(i,j,k,ibtcr)=                                          &
     &          MIN( MAX( bed(i,j,k,ibtcr), tcr_min(ng)), tcr_max(ng) )
          ENDDO
          ENDDO I_LOOP_CB2
#endif /* cohesive bed */
      END DO J_LOOP2
!
!-----------------------------------------------------------------------
! Store old bed thickness.
!-----------------------------------------------------------------------
!
# if defined SED_MORPH
      DO j=JstrR,JendR
        DO i=IstrR,IendR
            bed_thick(i,j,nnew)=0.0_r8
            DO kbed=1,Nbed
              bed_thick(i,j,nnew)=bed_thick(i,j,nnew)+                  &
     &                            bed(i,j,kbed,ithck)
            END DO
          END DO
        END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bed_thick(:,:,nnew))
#  endif
# endif
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
     &                    NghostPoints, EWperiodic, NSperiodic,         &
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
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed)
# endif

      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, MBOTP,                   &
     &                  bottom)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bottom)
# endif

      RETURN
      END SUBROUTINE bed_cohesive_tile
#endif
      END MODULE bed_cohesive_mod
