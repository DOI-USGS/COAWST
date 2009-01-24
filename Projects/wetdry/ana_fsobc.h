      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id: ana_fsobc.h 168 2008-03-18 20:12:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff, fac, omega, phase, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined WETDRY_SLOPE_CHAN
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=10.0_r8*sin(2.0_r8*pi/1.0_r8*      &
     &                               tdays(ng))-10.0_r8+0.10_r8
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
