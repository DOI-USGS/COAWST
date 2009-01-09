#include "cppdefs.h"

      MODULE bottom_mod
#if defined NONLINEAR && defined SEDIMENT
!
!svn $Id: bottom.F 1618 2008-08-03 03:05:47Z jcwarner $
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
      PUBLIC  :: bottom

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bottom (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
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
      CALL bottom_tile (ng, tile,                                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng), nnew(ng),                           &
     &                    OCEAN(ng) % bed_frac,                         &
     &                    OCEAN(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE bottom
!
!***********************************************************************
      SUBROUTINE bottom_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp, nnew,                             &
     &                          bed_frac, bottom)
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
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
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
      integer :: i, ised, j

      real(r8) :: cff1, cff2, cff3, cff4
      real(r8), parameter :: eps = 1.0E-14_r8

# include "set_bounds.h"

!
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines
!
      J_LOOP : DO j=Jstr,Jend
        DO i=Istr,Iend
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised,ng)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=MIN(cff2,Zob(ng))
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=MAX(cff4,1050.0_r8)
        END DO
      END DO J_LOOP

!
!  Move bottom roughness computations from BBL routines to here.
!

!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
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
      END SUBROUTINE bottom_tile
#endif
      END MODULE bottom_mod
