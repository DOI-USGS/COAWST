#include "cppdefs.h"
      MODULE exchange_2d_mod
!
!svn $Id: exchange_2d.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  2D fields.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                      Nested grid number.                      !
!     tile                    Domain partition.                        !
!     LBi                     I-dimension Lower bound.                 !
!     UBi                     I-dimension Upper bound.                 !
!     LBj                     J-dimension Lower bound.                 !
!     UBj                     J-dimension Upper bound.                 !
!     A                       2D field.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                       Processed 2D field                       !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     exchange_p2d_tile       periodic conditions at PSI-points        !
!     exchange_r2d_tile       periodic conditions at RHO-points        !
!     exchange_u2d_tile       periodic conditions at U-points          !
!     exchange_v2d_tile       periodic conditions at V-points          !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
#ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
#else
        EW_exchange=.TRUE.
#endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
#ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
#else
        NS_exchange=.TRUE.
#endif
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
            DO j=Jmin,Jmax
              A(Lm(ng)+1,j)=A(1,j)
              A(Lm(ng)+2,j)=A(2,j)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO j=Jmin,Jmax
                A(Lm(ng)+3,j)=A(3,j)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jmin,Jmax
              A(-2,j)=A(Lm(ng)-2,j)
              A(-1,j)=A(Lm(ng)-1,j)
              A( 0,j)=A(Lm(ng)  ,j)
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
            DO i=Imin,Imax
              A(i,Mm(ng)+1)=A(i,1)
              A(i,Mm(ng)+2)=A(i,2)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO i=Imin,Imax
                A(i,Mm(ng)+3)=A(i,3)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Imin,Imax
              A(i,-2)=A(i,Mm(ng)-2)
              A(i,-1)=A(i,Mm(ng)-1)
              A(i, 0)=A(i,Mm(ng)  )
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
            A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
            A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
            A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
            A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
              A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
              A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
              A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
              A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
            A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
            A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
            A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
            A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
            A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
            IF (NghostPoints.eq.3) THEN
              A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
              A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
              A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
            A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
            A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
            A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
            A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
            A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
              A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
              A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
            A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
            A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
            A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
            A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
            A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
            A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
            A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
            A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
#ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
#else
        EW_exchange=.TRUE.
#endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
#ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
#else
        NS_exchange=.TRUE.
#endif
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
            DO j=Jmin,Jmax
              A(Lm(ng)+1,j)=A(1,j)
              A(Lm(ng)+2,j)=A(2,j)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO j=Jmin,Jmax
                A(Lm(ng)+3,j)=A(3,j)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jmin,Jmax
              A(-2,j)=A(Lm(ng)-2,j)
              A(-1,j)=A(Lm(ng)-1,j)
              A( 0,j)=A(Lm(ng)  ,j)
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
            DO i=Imin,Imax
              A(i,Mm(ng)+1)=A(i,1)
              A(i,Mm(ng)+2)=A(i,2)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO i=Imin,Imax
                A(i,Mm(ng)+3)=A(i,3)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Imin,Imax
              A(i,-2)=A(i,Mm(ng)-2)
              A(i,-1)=A(i,Mm(ng)-1)
              A(i, 0)=A(i,Mm(ng)  )
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
            A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
            A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
            A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
            A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
              A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
              A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
              A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
              A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
            A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
            A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
            A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
            A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
            A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
            IF (NghostPoints.eq.3) THEN
              A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
              A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
              A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
            A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
            A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
            A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
            A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
            A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
              A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
              A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
            A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
            A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
            A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
            A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
            A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
            A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
            A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
            A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
#ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
#else
        EW_exchange=.TRUE.
#endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
#ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
#else
        NS_exchange=.TRUE.
#endif
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
            DO j=Jmin,Jmax
              A(Lm(ng)+1,j)=A(1,j)
              A(Lm(ng)+2,j)=A(2,j)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO j=Jmin,Jmax
                A(Lm(ng)+3,j)=A(3,j)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jmin,Jmax
              A(-2,j)=A(Lm(ng)-2,j)
              A(-1,j)=A(Lm(ng)-1,j)
              A( 0,j)=A(Lm(ng)  ,j)
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
            DO i=Imin,Imax
              A(i,Mm(ng)+1)=A(i,1)
              A(i,Mm(ng)+2)=A(i,2)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO i=Imin,Imax
                A(i,Mm(ng)+3)=A(i,3)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Imin,Imax
              A(i,-2)=A(i,Mm(ng)-2)
              A(i,-1)=A(i,Mm(ng)-1)
              A(i, 0)=A(i,Mm(ng)  )
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
            A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
            A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
            A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
            A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
              A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
              A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
              A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
              A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
            A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
            A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
            A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
            A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
            A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
            IF (NghostPoints.eq.3) THEN
              A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
              A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
              A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
            A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
            A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
            A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
            A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
            A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
              A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
              A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
            A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
            A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
            A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
            A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
            A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
            A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
            A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
            A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
#ifdef DISTRIBUTE
        EW_exchange=NtileI(ng).eq.1
#else
        EW_exchange=.TRUE.
#endif
      ELSE
        EW_exchange=.FALSE.
      END IF

      IF (NSperiodic(ng)) THEN
#ifdef DISTRIBUTE
        NS_exchange=NtileJ(ng).eq.1
#else
        NS_exchange=.TRUE.
#endif
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
            DO j=Jmin,Jmax
              A(Lm(ng)+1,j)=A(1,j)
              A(Lm(ng)+2,j)=A(2,j)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO j=Jmin,Jmax
                A(Lm(ng)+3,j)=A(3,j)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jmin,Jmax
              A(-2,j)=A(Lm(ng)-2,j)
              A(-1,j)=A(Lm(ng)-1,j)
              A( 0,j)=A(Lm(ng)  ,j)
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
            DO i=Imin,Imax
              A(i,Mm(ng)+1)=A(i,1)
              A(i,Mm(ng)+2)=A(i,2)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO i=Imin,Imax
                A(i,Mm(ng)+3)=A(i,3)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Imin,Imax
              A(i,-2)=A(i,Mm(ng)-2)
              A(i,-1)=A(i,Mm(ng)-1)
              A(i, 0)=A(i,Mm(ng)  )
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
            A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
            A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
            A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
            A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
              A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
              A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
              A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
              A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
            A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
            A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
            A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
            A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
            A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
            IF (NghostPoints.eq.3) THEN
              A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
              A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
              A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
            A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
            A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
            A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
            A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
            A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
            IF (NghostPoints.eq.3) THEN
              A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
              A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
              A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
            A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
            A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
            A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
            A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
            A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
            A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
            A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
            A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_v2d_tile

      END MODULE exchange_2d_mod
