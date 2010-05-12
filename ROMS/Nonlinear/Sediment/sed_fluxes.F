#include "cppdefs.h"

      MODULE sed_fluxes_mod

#if defined NONLINEAR && defined SEDIMENT && defined SUSPLOAD
!
!svn $Id: sed_fluxes.F 429 2009-12-20 17:30:26Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This computes sediment bed and water column exchanges: deposition,  !
!  resuspension, and erosion.                                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_fluxes

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_fluxes (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
      USE mod_bbl
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
      CALL sed_fluxes_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
     &                      GRID(ng) % Hz,                              &
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef BBL_MODEL
     &                      BBL(ng) % bustrc,                           &
     &                      BBL(ng) % bvstrc,                           &
     &                      BBL(ng) % bustrw,                           &
     &                      BBL(ng) % bvstrw,                           &
     &                      BBL(ng) % bustrcwmax,                       &
     &                      BBL(ng) % bvstrcwmax,                       &
# endif
     &                      FORCES(ng) % bustr,                         &
     &                      FORCES(ng) % bvstr,                         &
     &                      OCEAN(ng) % t,                              &
     &                      SEDBED(ng) % ero_flux,                      &
     &                      SEDBED(ng) % settling_flux,                 &
# if defined SED_MORPH
     &                      SEDBED(ng) % bed_thick,                     &
# endif
     &                      SEDBED(ng) % bed,                           &
     &                      SEDBED(ng) % bed_frac,                      &
     &                      SEDBED(ng) % bed_mass,                      &
     &                      SEDBED(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_fluxes
!
!***********************************************************************
      SUBROUTINE sed_fluxes_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
     &                            Hz,                                   &
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
# ifdef BBL_MODEL
     &                            bustrc, bvstrc,                       &
     &                            bustrw, bvstrw,                       &
     &                            bustrcwmax, bvstrcwmax,               &
# endif
     &                            bustr, bvstr,                         &
     &                            t,                                    &
     &                            ero_flux, settling_flux,              &
# if defined SED_MORPH
     &                            bed_thick,                            &
# endif
     &                            bed, bed_frac, bed_mass,              &
     &                            bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
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
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ero_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: settling_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
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
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ero_flux(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: settling_flux(LBi:UBi,LBj:UBj,NST)
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

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_w

# include "set_bounds.h"

# ifdef BEDLOAD
      bnew=nnew
# else
      bnew=nstp
# endif
!
!-----------------------------------------------------------------------
!  Compute sediment deposition, resuspension, and erosion.
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
!-----------------------------------------------------------------------
!  Sediment deposition and resuspension near the bottom.
!-----------------------------------------------------------------------
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precepitation settling_flux, already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
!
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
          DO i=Istr,Iend
!
!  Calculate critical shear stress in Pa
!
# if defined COHESIVE_BED
            cff = rho0/bed(i,j,1,ibtcr)
# elif defined MIXED_BED
            cff = MAX(bottom(i,j,idprp)*bed(i,j,1,ibtcr)/rho0+          &
     &            (1.0_r8-bottom(i,j,idprp))*tau_ce(ised,ng),           &
     &            tau_ce(ised,ng))
            cff=1.0_r8/cff
# else
            cff=1.0_r8/tau_ce(ised,ng)
# endif
!
!  Compute erosion, ero_flux (kg/m2).
!
            cff1=(1.0_r8-bed(i,j,1,iporo))*bed_frac(i,j,1,ised)
            cff2=dt(ng)*Erate(ised,ng)*cff1
            cff3=Srho(ised,ng)*cff1
            cff4=bed_mass(i,j,1,bnew,ised)
            ero_flux(i,j,ised)=                                         &
     &               MIN(MAX(0.0_r8,cff2*(cff*tau_w(i,j)-1.0_r8)),      &
     &                   MIN(cff3*bottom(i,j,iactv),cff4)+              &
     &                   settling_flux(i,j,ised))
!
!  Update global tracer variables (m Tunits for nnew indx, Tuints for 3)
!  for erosive flux.
!
            t(i,j,1,nnew,indx)=t(i,j,1,nnew,indx)+ero_flux(i,j,ised)
# ifdef TS_MPDATA
            t(i,j,1,3,indx)=t(i,j,1,3,indx)+                            &
     &                      ero_flux(i,j,ised)*Hz_inv(i,1)
# endif
          END DO
        END DO SED_LOOP
      END DO J_LOOP

      RETURN
      END SUBROUTINE sed_fluxes_tile
#endif
      END MODULE sed_fluxes_mod
