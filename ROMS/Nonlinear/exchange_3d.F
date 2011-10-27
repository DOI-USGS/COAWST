#include "cppdefs.h"
      MODULE exchange_3d_mod
#if defined SOLVE3D && (defined EW_PERIODIC || defined NS_PERIODIC)
!
!svn $Id: exchange_3d.F 1597 2008-07-25 02:31:35Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This package contains periodic boundary conditions and parallel     !
!  exchage (distributed-memory only) routines for 3D variables.        !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p3d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r3d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u3d_tile    periodic conditions/exchange at U-points     !
!    exchange_v3d_tile    periodic conditions/exchange at V-points     !
!    exchange_w3d_tile    periodic conditions/exchange at W-points     !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
# ifdef REFINED_GRID
!  Set this for Parent Grid only.
      IF (ng.eq.1) THEN
# endif
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
#  ifdef THREE_GHOST
              A(Lm(ng)+3,j,k)=A(3,j,k)
#  endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
#  ifdef THREE_GHOST
              A(i,Mm(ng)+3,k)=A(i,3,k)
#  endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
#  endif
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
            A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
            A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
            A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
#  ifdef THREE_GHOST
            A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
            A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
            A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
#  endif
          END DO
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
            A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
            A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_p3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
# ifdef REFINED_GRID
!  Set this for Parent Grid only.
      IF (ng.eq.1) THEN
# endif
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
#  ifdef THREE_GHOST
              A(Lm(ng)+3,j,k)=A(3,j,k)
#  endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
#  ifdef THREE_GHOST
              A(i,Mm(ng)+3,k)=A(i,3,k)
#  endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif
# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
#  endif
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
            A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
            A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
            A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
#  ifdef THREE_GHOST
            A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
            A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
            A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
#  endif
          END DO
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
            A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
            A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_r3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
# ifdef REFINED_GRID
!  Set this for Parent Grid only.
      IF (ng.eq.1) THEN
# endif
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
#  ifdef THREE_GHOST
              A(Lm(ng)+3,j,k)=A(3,j,k)
#  endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
#  ifdef THREE_GHOST
              A(i,Mm(ng)+3,k)=A(i,3,k)
#  endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
#  endif
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
            A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
            A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
            A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
#  ifdef THREE_GHOST
            A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
            A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
            A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
#  endif
          END DO
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
            A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
            A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_u3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
# ifdef REFINED_GRID
!  Set this for Parent Grid only.
      IF (ng.eq.1) THEN
# endif
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
#  ifdef THREE_GHOST
              A(Lm(ng)+3,j,k)=A(3,j,k)
#  endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
#  ifdef THREE_GHOST
              A(i,Mm(ng)+3,k)=A(i,3,k)
#  endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
#  endif
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
            A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
            A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
            A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
#  ifdef THREE_GHOST
            A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
            A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
            A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
#  endif
          END DO
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
            A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
            A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_v3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
# ifdef REFINED_GRID
!  Set this for Parent Grid only.
      IF (ng.eq.1) THEN
# endif
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
#  ifdef THREE_GHOST
              A(Lm(ng)+3,j,k)=A(3,j,k)
#  endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
#  ifdef THREE_GHOST
              A(i,Mm(ng)+3,k)=A(i,3,k)
#  endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif
# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
#  endif
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
            A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
            A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
            A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
#  ifdef THREE_GHOST
            A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
            A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
            A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
#  endif
          END DO
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
            A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
            A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
#  endif
          END DO
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif
# ifdef REFINED_GRID
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_w3d_tile
#endif

      END MODULE exchange_3d_mod
