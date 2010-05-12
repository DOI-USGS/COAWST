      SUBROUTINE ana_m3obc (ng, tile, model)
!
!! svn $Id: ana_m3obc.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
      CALL ana_m3obc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(14)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m3obc
!
!***********************************************************************
      SUBROUTINE ana_m3obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  3D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
# ifdef EAST_M3OBC
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
# endif
# ifdef WEST_M3OBC
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
# endif
# ifdef SOUTH_M3OBC
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
# endif
# ifdef NORTH_M3OBC
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
# endif
#else
      ana_m3obc.h: No values provided for BOUNDARY(ng)%u_xxxx and
                                          BOUNDARY(ng)%v_xxxx
#endif

      RETURN
      END SUBROUTINE ana_m3obc_tile
