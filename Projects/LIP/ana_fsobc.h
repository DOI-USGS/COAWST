      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id: ana_fsobc.h 737 2008-09-07 02:06:44Z jcwarner $
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
#if defined SHOREFACE
      IF (WESTERN_EDGE) THEN
        cff=-4.1_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
#else
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
