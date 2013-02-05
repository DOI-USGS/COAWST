#include "cppdefs.h"

      MODULE sed_biomass_mod

#if defined NONLINEAR && defined SEDIMENT && defined SED_BIOMASS
!
!svn $Id: sed_biomass.F 429 2009-12-20 17:30:26Z jcwarner $
!==================================================== Neil K. Ganju ====
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This computes sediment biomass due to vegation growth               !
!                                                                      !
!  References:                                                         !
!                                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_biomass

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_biomass (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sed_biomass_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
#  endif
# endif
     &                      GRID(ng) % z_w,                             &
     &                      SEDBED(ng) % Dstp_max,                      &
     &                      SEDBED(ng) % bottom,                        &
     &                      SEDBED(ng) % settling_flux)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_biomass
!
!***********************************************************************
      SUBROUTINE sed_biomass_tile (ng, tile,                            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
# ifdef MASKING
     &                            rmask,                                &
#  ifdef WET_DRY
     &                            rmask_wet,                            &
#  endif
# endif
     &                            z_w,                                  &
     &                            Dstp_max,                             &
     &                            bottom, settling_flux)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Dstp_max(LBi:,LBj:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
      real(r8), intent(inout) :: settling_flux(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: Dstp_max(LBi:UBi,LBj:UBj,MBOTP)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
      real(r8), intent(inout) :: settling_flux(LBi:UBi,LBj:UBj,NST)
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
      integer :: i, j, k, ised
      integer :: sstp, nbio_steps

      real(r8) :: cff, Dstp

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Compute 
!-----------------------------------------------------------------------
!
!  Compute number of model steps for each hour.
!
      nbio_steps=MAX(1,INT(3600.0_r8/dt(ng)))
!
!  Compute number of hourly values we need to save.
!  If we want a 1 day avg, then need 24 values.
!
      J_LOOP : DO j=Jstr,Jend
!
!  Only update the max depth once per hour.
!
        IF (MOD(iic(ng),nbio_steps).eq.0) THEN
          DO i=Istr,Iend
!
!  Determine the index for placement of new value.
!
            sstp=1+MOD(iic(ng)-ntstart(ng),24)
!
!  Save instantaneous depth at this instance and recompute max daily depth.
!
            Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
            Dstp_max(i,j,sstp)=Dstp
            cff=0.0_r8
            DO k=1,24                        ! should be a variable, also above
              cff=MAX(cff,Dstp_max(i,j,k))
            END DO
            bottom(i,j,imaxD)=cff
          END DO
        END IF
!
!  Update settling flux for depositing bio mass.
!
        DO i=Istr,Iend
!
!  Require (for now) that the first sed class be the new biomass.
!
          ised=1
          cff=0.0_r8                               ! remove this line
!         cff= funct( bottom(i,j,imaxD), dt(ng))   ! need real eq. in kg/m^2
# ifdef MASKING
          cff=cff*rmask(i,j)
#  ifdef WET_DRY
          cff=cff*rmask_wet(i,j)          ! not sure this is needed
#  endif
# endif
          settling_flux(i,j,ised)=settling_flux(i,j,ised)+cff
        END DO
      END DO J_LOOP

      RETURN
      END SUBROUTINE sed_biomass_tile
#endif
      END MODULE sed_biomass_mod
