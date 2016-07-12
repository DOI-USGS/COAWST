#include "cppdefs.h"
      MODULE exchange_3d_mod
#ifdef SOLVE3D
!
!svn $Id: exchange_3d.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  3D fields.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                      Nested grid number.                      !
!     tile                    Domain partition.                        !
!     LBi                     I-dimension Lower bound.                 !
!     UBi                     I-dimension Upper bound.                 !
!     LBj                     J-dimension Lower bound.                 !
!     UBj                     J-dimension Upper bound.                 !
!     LBk                     J-dimension Lower bound.                 !
!     UBk                     J-dimension Upper bound.                 !
!     A                       3D field.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                       Processed 3D field.                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     exchange_p3d_tile       periodic conditions at PSI-points        !
!     exchange_r3d_tile       periodic conditions at RHO-points        !
!     exchange_u3d_tile       periodic conditions at U-points          !
!     exchange_v3d_tile       periodic conditions at V-points          !
!     exchange_w3d_tile       periodic conditions at W-points          !
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
      USE mod_scalars
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
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=Jstr
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=Istr
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

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
      USE mod_scalars
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
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

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
      USE mod_scalars
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
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=Istr
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

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
      USE mod_scalars
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
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=Jstr
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

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
      USE mod_scalars
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
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_w3d_tile

# ifdef INWAVE_MODEL
!
!***********************************************************************
      SUBROUTINE exchange_AC3d_tile (ng, tile,                          &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
#  else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
# ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
# else
        EW_exchange=.TRUE.
# endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
# ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
# else
        NS_exchange=.TRUE.
# endif
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_AC3d_tile
# endif

#endif

      END MODULE exchange_3d_mod
