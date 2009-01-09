      SUBROUTINE ana_bmflux (ng, tile, model)
!
!! svn $Id: ana_bmflux.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets spatially varying bottom roughness Zo, rdrg2, or  !
!  rdrg parameter using an analytical expression.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_ocean
      USE mod_grid
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_bmflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
     &                      OCEAN(ng) % bottom)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME( 2)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_bmflux
!
!***********************************************************************
      SUBROUTINE ana_bmflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            xr, yr,                               &
     &                            bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment

#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
      real(r8), intent(out) :: bottom(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set spatially varying bottom Zo.
!-----------------------------------------------------------------------
!
#if defined SED_TOY
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,izdef)=0.002_r8
        END DO
      END DO
#else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,izdef)=Zob(ng)
        END DO
      END DO
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, MBOTP,             &
     &                        NghostPoints,                             &
     &                        bottom)
#endif
      RETURN
      END SUBROUTINE ana_bmflux_tile
