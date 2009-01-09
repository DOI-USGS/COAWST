      SUBROUTINE ana_m3obc (ng, tile, model)
!
!! svn $Id: ana_m3obc.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets 3D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m3obc_tile (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(14),'(a,a)') TRIM(Adir), '/ana_m3obc.h'
      END IF

      RETURN
      END SUBROUTINE ana_m3obc
!
!***********************************************************************
      SUBROUTINE ana_m3obc_tile (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  3D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_east(j,k)=???
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_east(j,k)=???
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_west(j,k)=???
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_west(j,k)=???
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_south(i,k)=???
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_south(i,k)=???
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_north(i,k)=???
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_north(i,k)=???
          END DO
        END DO
      END IF
#else
      ana_m3obc.h: No values provided for BOUNDARY(ng)%u_xxxx and
                                          BOUNDARY(ng)%v_xxxx
#endif

      RETURN
      END SUBROUTINE ana_m3obc_tile
