      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id: ana_fsobc.h 735 2007-04-27 14:00:46Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
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
      CALL ana_fsobc_tile (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME( 6),'(a,a)') TRIM(Adir), '/ana_fsobc.h'
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, model, Istr, Iend, Jstr, Jend,     &
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
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
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
#elif defined VISSER
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF
#else
      ana_fsobc_user.h: No values provided for BOUNDARY(ng)%zeta_xxxx.
#endif

      RETURN
      END SUBROUTINE ana_fsobc_tile
