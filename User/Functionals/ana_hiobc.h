      SUBROUTINE ana_hiobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
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
      CALL ana_hiobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(44)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_hiobc
!
!***********************************************************************
      SUBROUTINE ana_hiobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
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
      IF (LBC(iwest,isHice,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%hi_west(j)=0.0_r8
        END DO
      END IF
      IF (LBC(ieast,isHice,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%hi_east(j)=0.0_r8
        END DO
      END IF
      IF (LBC(isouth,isHice,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%hi_south(i)=0.0_r8
        END DO
      END IF
      IF (LBC(inorth,isHice,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%hi_north(i)=0.0_r8
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_hiobc_tile
