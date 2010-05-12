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
      real(r8) :: fac, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  3D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined SED_TEST1
# ifdef WEST_M3OBC
      IF (WESTERN_EDGE) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrR,JendR
            val=0.5_r8*(zeta(0 ,j,knew)+h(0 ,j)+                        &
     &                  zeta(1 ,j,knew)+h(1 ,j))
            BOUNDARY(ng)%u_west(j,k)=-LOG((val+0.5*(z_r(Istr-1,j,k)+    &
     &                                              z_r(Istr  ,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef EAST_M3OBC
      IF (EASTERN_EDGE) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrR,JendR
            val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                &
     &                  zeta(Iend+1,j,knew)+h(Iend+1,j))
            BOUNDARY(ng)%u_east(j,k)=-LOG((val+0.5*(z_r(Iend  ,j,k)+    &
     &                                              z_r(Iend+1,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
#else
# ifdef EAST_M3OBC
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_east(j,k)=0.0_r8
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef WEST_M3OBC
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_west(j,k)=0.0_r8
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef SOUTH_M3OBC
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_south(i,k)=0.0_r8
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_south(i,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef NORTH_M3OBC
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_north(i,k)=0.0_r8
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_north(i,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
#endif
      RETURN
      END SUBROUTINE ana_m3obc_tile
