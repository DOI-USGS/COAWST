#include "cppdefs.h"
      MODULE exchange_4d_mod
#ifdef SOLVE3D
!
!svn $Id: exchange_4d.F 1054 2021-03-06 19:47:12Z arango $
!=======================================================================
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!     LBl                     L-dimension Lower bound.                 !
!     UBl                     L-dimension Upper bound.                 !
!     A                       4D field.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                       Processed 3D field.                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     exchange_p4d_tile       periodic conditions at PSI-points        !
!     exchange_r4d_tile       periodic conditions at RHO-points        !
!     exchange_u4d_tile       periodic conditions at U-points          !
!     exchange_v4d_tile       periodic conditions at V-points          !
!     exchange_w4d_tile       periodic conditions at W-points          !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p4d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              LBl, UBl,                           &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+1,j,k,l)=A(1,j,k,l)
                  A(Lm(ng)+2,j,k,l)=A(2,j,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO j=Jmin,Jmax
                    A(Lm(ng)+3,j,k,l)=A(3,j,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(-2,j,k,l)=A(Lm(ng)-2,j,k,l)
                  A(-1,j,k,l)=A(Lm(ng)-1,j,k,l)
                  A( 0,j,k,l)=A(Lm(ng)  ,j,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+1,k,l)=A(i,1,k,l)
                  A(i,Mm(ng)+2,k,l)=A(i,2,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO i=Imin,Imax
                    A(i,Mm(ng)+3,k,l)=A(i,3,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,-2,k,l)=A(i,Mm(ng)-2,k,l)
                  A(i,-1,k,l)=A(i,Mm(ng)-1,k,l)
                  A(i, 0,k,l)=A(i,Mm(ng)  ,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+1,k,l)=A(1,1,k,l)
                A(Lm(ng)+1,Mm(ng)+2,k,l)=A(1,2,k,l)
                A(Lm(ng)+2,Mm(ng)+1,k,l)=A(2,1,k,l)
                A(Lm(ng)+2,Mm(ng)+2,k,l)=A(2,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+1,Mm(ng)+3,k,l)=A(1,3,k,l)
                  A(Lm(ng)+2,Mm(ng)+3,k,l)=A(2,3,k,l)
                  A(Lm(ng)+3,Mm(ng)+1,k,l)=A(3,1,k,l)
                  A(Lm(ng)+3,Mm(ng)+2,k,l)=A(3,2,k,l)
                  A(Lm(ng)+3,Mm(ng)+3,k,l)=A(3,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,Mm(ng)+1,k,l)=A(Lm(ng)-2,1,k,l)
                A(-1,Mm(ng)+1,k,l)=A(Lm(ng)-1,1,k,l)
                A( 0,Mm(ng)+1,k,l)=A(Lm(ng)  ,1,k,l)
                A(-2,Mm(ng)+2,k,l)=A(Lm(ng)-2,2,k,l)
                A(-1,Mm(ng)+2,k,l)=A(Lm(ng)-1,2,k,l)
                A( 0,Mm(ng)+2,k,l)=A(Lm(ng)  ,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(-2,Mm(ng)+3,k,l)=A(Lm(ng)-2,3,k,l)
                  A(-1,Mm(ng)+3,k,l)=A(Lm(ng)-1,3,k,l)
                  A( 0,Mm(ng)+3,k,l)=A(Lm(ng)  ,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,-2,k,l)=A(1,Mm(ng)-2,k,l)
                A(Lm(ng)+1,-1,k,l)=A(1,Mm(ng)-1,k,l)
                A(Lm(ng)+1, 0,k,l)=A(1,Mm(ng)  ,k,l)
                A(Lm(ng)+2,-2,k,l)=A(2,Mm(ng)-2,k,l)
                A(Lm(ng)+2,-1,k,l)=A(2,Mm(ng)-1,k,l)
                A(Lm(ng)+2, 0,k,l)=A(2,Mm(ng)  ,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+3,-2,k,l)=A(3,Mm(ng)-2,k,l)
                  A(Lm(ng)+3,-1,k,l)=A(3,Mm(ng)-1,k,l)
                  A(Lm(ng)+3, 0,k,l)=A(3,Mm(ng)  ,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,-2,k,l)=A(Lm(ng)-2,Mm(ng)-2,k,l)
                A(-2,-1,k,l)=A(Lm(ng)-2,Mm(ng)-1,k,l)
                A(-2, 0,k,l)=A(Lm(ng)-2,Mm(ng)  ,k,l)
                A(-1,-2,k,l)=A(Lm(ng)-1,Mm(ng)-2,k,l)
                A(-1,-1,k,l)=A(Lm(ng)-1,Mm(ng)-1,k,l)
                A(-1, 0,k,l)=A(Lm(ng)-1,Mm(ng)  ,k,l)
                A( 0,-2,k,l)=A(Lm(ng)  ,Mm(ng)-2,k,l)
                A( 0,-1,k,l)=A(Lm(ng)  ,Mm(ng)-1,k,l)
                A( 0, 0,k,l)=A(Lm(ng)  ,Mm(ng)  ,k,l)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_p4d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r4d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              LBl, UBl,                           &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+1,j,k,l)=A(1,j,k,l)
                  A(Lm(ng)+2,j,k,l)=A(2,j,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO j=Jmin,Jmax
                    A(Lm(ng)+3,j,k,l)=A(3,j,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(-2,j,k,l)=A(Lm(ng)-2,j,k,l)
                  A(-1,j,k,l)=A(Lm(ng)-1,j,k,l)
                  A( 0,j,k,l)=A(Lm(ng)  ,j,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+1,k,l)=A(i,1,k,l)
                  A(i,Mm(ng)+2,k,l)=A(i,2,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO i=Imin,Imax
                    A(i,Mm(ng)+3,k,l)=A(i,3,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,-2,k,l)=A(i,Mm(ng)-2,k,l)
                  A(i,-1,k,l)=A(i,Mm(ng)-1,k,l)
                  A(i, 0,k,l)=A(i,Mm(ng)  ,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+1,k,l)=A(1,1,k,l)
                A(Lm(ng)+1,Mm(ng)+2,k,l)=A(1,2,k,l)
                A(Lm(ng)+2,Mm(ng)+1,k,l)=A(2,1,k,l)
                A(Lm(ng)+2,Mm(ng)+2,k,l)=A(2,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+1,Mm(ng)+3,k,l)=A(1,3,k,l)
                  A(Lm(ng)+2,Mm(ng)+3,k,l)=A(2,3,k,l)
                  A(Lm(ng)+3,Mm(ng)+1,k,l)=A(3,1,k,l)
                  A(Lm(ng)+3,Mm(ng)+2,k,l)=A(3,2,k,l)
                  A(Lm(ng)+3,Mm(ng)+3,k,l)=A(3,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,Mm(ng)+1,k,l)=A(Lm(ng)-2,1,k,l)
                A(-1,Mm(ng)+1,k,l)=A(Lm(ng)-1,1,k,l)
                A( 0,Mm(ng)+1,k,l)=A(Lm(ng)  ,1,k,l)
                A(-2,Mm(ng)+2,k,l)=A(Lm(ng)-2,2,k,l)
                A(-1,Mm(ng)+2,k,l)=A(Lm(ng)-1,2,k,l)
                A( 0,Mm(ng)+2,k,l)=A(Lm(ng)  ,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(-2,Mm(ng)+3,k,l)=A(Lm(ng)-2,3,k,l)
                  A(-1,Mm(ng)+3,k,l)=A(Lm(ng)-1,3,k,l)
                  A( 0,Mm(ng)+3,k,l)=A(Lm(ng)  ,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,-2,k,l)=A(1,Mm(ng)-2,k,l)
                A(Lm(ng)+1,-1,k,l)=A(1,Mm(ng)-1,k,l)
                A(Lm(ng)+1, 0,k,l)=A(1,Mm(ng)  ,k,l)
                A(Lm(ng)+2,-2,k,l)=A(2,Mm(ng)-2,k,l)
                A(Lm(ng)+2,-1,k,l)=A(2,Mm(ng)-1,k,l)
                A(Lm(ng)+2, 0,k,l)=A(2,Mm(ng)  ,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+3,-2,k,l)=A(3,Mm(ng)-2,k,l)
                  A(Lm(ng)+3,-1,k,l)=A(3,Mm(ng)-1,k,l)
                  A(Lm(ng)+3, 0,k,l)=A(3,Mm(ng)  ,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,-2,k,l)=A(Lm(ng)-2,Mm(ng)-2,k,l)
                A(-2,-1,k,l)=A(Lm(ng)-2,Mm(ng)-1,k,l)
                A(-2, 0,k,l)=A(Lm(ng)-2,Mm(ng)  ,k,l)
                A(-1,-2,k,l)=A(Lm(ng)-1,Mm(ng)-2,k,l)
                A(-1,-1,k,l)=A(Lm(ng)-1,Mm(ng)-1,k,l)
                A(-1, 0,k,l)=A(Lm(ng)-1,Mm(ng)  ,k,l)
                A( 0,-2,k,l)=A(Lm(ng)  ,Mm(ng)-2,k,l)
                A( 0,-1,k,l)=A(Lm(ng)  ,Mm(ng)-1,k,l)
                A( 0, 0,k,l)=A(Lm(ng)  ,Mm(ng)  ,k,l)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_r4d_tile
!
!***********************************************************************
      SUBROUTINE exchange_u4d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              LBl, UBl,                           &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+1,j,k,l)=A(1,j,k,l)
                  A(Lm(ng)+2,j,k,l)=A(2,j,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO j=Jmin,Jmax
                    A(Lm(ng)+3,j,k,l)=A(3,j,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(-2,j,k,l)=A(Lm(ng)-2,j,k,l)
                  A(-1,j,k,l)=A(Lm(ng)-1,j,k,l)
                  A( 0,j,k,l)=A(Lm(ng)  ,j,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+1,k,l)=A(i,1,k,l)
                  A(i,Mm(ng)+2,k,l)=A(i,2,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO i=Imin,Imax
                    A(i,Mm(ng)+3,k,l)=A(i,3,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,-2,k,l)=A(i,Mm(ng)-2,k,l)
                  A(i,-1,k,l)=A(i,Mm(ng)-1,k,l)
                  A(i, 0,k,l)=A(i,Mm(ng)  ,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+1,k,l)=A(1,1,k,l)
                A(Lm(ng)+1,Mm(ng)+2,k,l)=A(1,2,k,l)
                A(Lm(ng)+2,Mm(ng)+1,k,l)=A(2,1,k,l)
                A(Lm(ng)+2,Mm(ng)+2,k,l)=A(2,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+2,Mm(ng)+3,k,l)=A(2,3,k,l)
                  A(Lm(ng)+3,Mm(ng)+1,k,l)=A(3,1,k,l)
                  A(Lm(ng)+3,Mm(ng)+2,k,l)=A(3,2,k,l)
                  A(Lm(ng)+3,Mm(ng)+3,k,l)=A(3,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,Mm(ng)+1,k,l)=A(Lm(ng)-2,1,k,l)
                A(-1,Mm(ng)+1,k,l)=A(Lm(ng)-1,1,k,l)
                A( 0,Mm(ng)+1,k,l)=A(Lm(ng)  ,1,k,l)
                A(-2,Mm(ng)+2,k,l)=A(Lm(ng)-2,2,k,l)
                A(-1,Mm(ng)+2,k,l)=A(Lm(ng)-1,2,k,l)
                A( 0,Mm(ng)+2,k,l)=A(Lm(ng)  ,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(-2,Mm(ng)+3,k,l)=A(Lm(ng)-2,3,k,l)
                  A(-1,Mm(ng)+3,k,l)=A(Lm(ng)-1,3,k,l)
                  A( 0,Mm(ng)+3,k,l)=A(Lm(ng)  ,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,-2,k,l)=A(1,Mm(ng)-2,k,l)
                A(Lm(ng)+1,-1,k,l)=A(1,Mm(ng)-1,k,l)
                A(Lm(ng)+1, 0,k,l)=A(1,Mm(ng)  ,k,l)
                A(Lm(ng)+2,-2,k,l)=A(2,Mm(ng)-2,k,l)
                A(Lm(ng)+2,-1,k,l)=A(2,Mm(ng)-1,k,l)
                A(Lm(ng)+2, 0,k,l)=A(2,Mm(ng)  ,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+3,-2,k,l)=A(3,Mm(ng)-2,k,l)
                  A(Lm(ng)+3,-1,k,l)=A(3,Mm(ng)-1,k,l)
                  A(Lm(ng)+3, 0,k,l)=A(3,Mm(ng)  ,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,-2,k,l)=A(Lm(ng)-2,Mm(ng)-2,k,l)
                A(-2,-1,k,l)=A(Lm(ng)-2,Mm(ng)-1,k,l)
                A(-2, 0,k,l)=A(Lm(ng)-2,Mm(ng)  ,k,l)
                A(-1,-2,k,l)=A(Lm(ng)-1,Mm(ng)-2,k,l)
                A(-1,-1,k,l)=A(Lm(ng)-1,Mm(ng)-1,k,l)
                A(-1, 0,k,l)=A(Lm(ng)-1,Mm(ng)  ,k,l)
                A( 0,-2,k,l)=A(Lm(ng)  ,Mm(ng)-2,k,l)
                A( 0,-1,k,l)=A(Lm(ng)  ,Mm(ng)-1,k,l)
                A( 0, 0,k,l)=A(Lm(ng)  ,Mm(ng)  ,k,l)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_u4d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v4d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              LBl, UBl,                           &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+1,j,k,l)=A(1,j,k,l)
                  A(Lm(ng)+2,j,k,l)=A(2,j,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO j=Jmin,Jmax
                    A(Lm(ng)+3,j,k,l)=A(3,j,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(-2,j,k,l)=A(Lm(ng)-2,j,k,l)
                  A(-1,j,k,l)=A(Lm(ng)-1,j,k,l)
                  A( 0,j,k,l)=A(Lm(ng)  ,j,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+1,k,l)=A(i,1,k,l)
                  A(i,Mm(ng)+2,k,l)=A(i,2,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO i=Imin,Imax
                    A(i,Mm(ng)+3,k,l)=A(i,3,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,-2,k,l)=A(i,Mm(ng)-2,k,l)
                  A(i,-1,k,l)=A(i,Mm(ng)-1,k,l)
                  A(i, 0,k,l)=A(i,Mm(ng)  ,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+1,k,l)=A(1,1,k,l)
                A(Lm(ng)+1,Mm(ng)+2,k,l)=A(1,2,k,l)
                A(Lm(ng)+2,Mm(ng)+1,k,l)=A(2,1,k,l)
                A(Lm(ng)+2,Mm(ng)+2,k,l)=A(2,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+1,Mm(ng)+3,k,l)=A(1,3,k,l)
                  A(Lm(ng)+2,Mm(ng)+3,k,l)=A(2,3,k,l)
                  A(Lm(ng)+3,Mm(ng)+1,k,l)=A(3,1,k,l)
                  A(Lm(ng)+3,Mm(ng)+2,k,l)=A(3,2,k,l)
                  A(Lm(ng)+3,Mm(ng)+3,k,l)=A(3,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,Mm(ng)+1,k,l)=A(Lm(ng)-2,1,k,l)
                A(-1,Mm(ng)+1,k,l)=A(Lm(ng)-1,1,k,l)
                A( 0,Mm(ng)+1,k,l)=A(Lm(ng)  ,1,k,l)
                A(-2,Mm(ng)+2,k,l)=A(Lm(ng)-2,2,k,l)
                A(-1,Mm(ng)+2,k,l)=A(Lm(ng)-1,2,k,l)
                A( 0,Mm(ng)+2,k,l)=A(Lm(ng)  ,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(-2,Mm(ng)+3,k,l)=A(Lm(ng)-2,3,k,l)
                  A(-1,Mm(ng)+3,k,l)=A(Lm(ng)-1,3,k,l)
                  A( 0,Mm(ng)+3,k,l)=A(Lm(ng)  ,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,-2,k,l)=A(1,Mm(ng)-2,k,l)
                A(Lm(ng)+1,-1,k,l)=A(1,Mm(ng)-1,k,l)
                A(Lm(ng)+1, 0,k,l)=A(1,Mm(ng)  ,k,l)
                A(Lm(ng)+2,-2,k,l)=A(2,Mm(ng)-2,k,l)
                A(Lm(ng)+2,-1,k,l)=A(2,Mm(ng)-1,k,l)
                A(Lm(ng)+2, 0,k,l)=A(2,Mm(ng)  ,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+3,-2,k,l)=A(3,Mm(ng)-2,k,l)
                  A(Lm(ng)+3,-1,k,l)=A(3,Mm(ng)-1,k,l)
                  A(Lm(ng)+3, 0,k,l)=A(3,Mm(ng)  ,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,-2,k,l)=A(Lm(ng)-2,Mm(ng)-2,k,l)
                A(-2,-1,k,l)=A(Lm(ng)-2,Mm(ng)-1,k,l)
                A(-2, 0,k,l)=A(Lm(ng)-2,Mm(ng)  ,k,l)
                A(-1,-2,k,l)=A(Lm(ng)-1,Mm(ng)-2,k,l)
                A(-1,-1,k,l)=A(Lm(ng)-1,Mm(ng)-1,k,l)
                A(-1, 0,k,l)=A(Lm(ng)-1,Mm(ng)  ,k,l)
                A( 0,-2,k,l)=A(Lm(ng)  ,Mm(ng)-2,k,l)
                A( 0,-1,k,l)=A(Lm(ng)  ,Mm(ng)-1,k,l)
                A( 0, 0,k,l)=A(Lm(ng)  ,Mm(ng)  ,k,l)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_v4d_tile

!
!***********************************************************************
      SUBROUTINE exchange_w4d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              LBl, UBl,                           &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+1,j,k,l)=A(1,j,k,l)
                  A(Lm(ng)+2,j,k,l)=A(2,j,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO j=Jmin,Jmax
                    A(Lm(ng)+3,j,k,l)=A(3,j,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(-2,j,k,l)=A(Lm(ng)-2,j,k,l)
                  A(-1,j,k,l)=A(Lm(ng)-1,j,k,l)
                  A( 0,j,k,l)=A(Lm(ng)  ,j,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+1,k,l)=A(i,1,k,l)
                  A(i,Mm(ng)+2,k,l)=A(i,2,k,l)
                END DO
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  DO i=Imin,Imax
                    A(i,Mm(ng)+3,k,l)=A(i,3,k,l)
                  END DO
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,-2,k,l)=A(i,Mm(ng)-2,k,l)
                  A(i,-1,k,l)=A(i,Mm(ng)-1,k,l)
                  A(i, 0,k,l)=A(i,Mm(ng)  ,k,l)
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+1,k,l)=A(1,1,k,l)
                A(Lm(ng)+1,Mm(ng)+2,k,l)=A(1,2,k,l)
                A(Lm(ng)+2,Mm(ng)+1,k,l)=A(2,1,k,l)
                A(Lm(ng)+2,Mm(ng)+2,k,l)=A(2,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+1,Mm(ng)+3,k,l)=A(1,3,k,l)
                  A(Lm(ng)+2,Mm(ng)+3,k,l)=A(2,3,k,l)
                  A(Lm(ng)+3,Mm(ng)+1,k,l)=A(3,1,k,l)
                  A(Lm(ng)+3,Mm(ng)+2,k,l)=A(3,2,k,l)
                  A(Lm(ng)+3,Mm(ng)+3,k,l)=A(3,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,Mm(ng)+1,k,l)=A(Lm(ng)-2,1,k,l)
                A(-1,Mm(ng)+1,k,l)=A(Lm(ng)-1,1,k,l)
                A( 0,Mm(ng)+1,k,l)=A(Lm(ng)  ,1,k,l)
                A(-2,Mm(ng)+2,k,l)=A(Lm(ng)-2,2,k,l)
                A(-1,Mm(ng)+2,k,l)=A(Lm(ng)-1,2,k,l)
                A( 0,Mm(ng)+2,k,l)=A(Lm(ng)  ,2,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(-2,Mm(ng)+3,k,l)=A(Lm(ng)-2,3,k,l)
                  A(-1,Mm(ng)+3,k,l)=A(Lm(ng)-1,3,k,l)
                  A( 0,Mm(ng)+3,k,l)=A(Lm(ng)  ,3,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Lm(ng)+1,-2,k,l)=A(1,Mm(ng)-2,k,l)
                A(Lm(ng)+1,-1,k,l)=A(1,Mm(ng)-1,k,l)
                A(Lm(ng)+1, 0,k,l)=A(1,Mm(ng)  ,k,l)
                A(Lm(ng)+2,-2,k,l)=A(2,Mm(ng)-2,k,l)
                A(Lm(ng)+2,-1,k,l)=A(2,Mm(ng)-1,k,l)
                A(Lm(ng)+2, 0,k,l)=A(2,Mm(ng)  ,k,l)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO l=LBl,UBl
                DO k=LBk,UBk
                  A(Lm(ng)+3,-2,k,l)=A(3,Mm(ng)-2,k,l)
                  A(Lm(ng)+3,-1,k,l)=A(3,Mm(ng)-1,k,l)
                  A(Lm(ng)+3, 0,k,l)=A(3,Mm(ng)  ,k,l)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(-2,-2,k,l)=A(Lm(ng)-2,Mm(ng)-2,k,l)
                A(-2,-1,k,l)=A(Lm(ng)-2,Mm(ng)-1,k,l)
                A(-2, 0,k,l)=A(Lm(ng)-2,Mm(ng)  ,k,l)
                A(-1,-2,k,l)=A(Lm(ng)-1,Mm(ng)-2,k,l)
                A(-1,-1,k,l)=A(Lm(ng)-1,Mm(ng)-1,k,l)
                A(-1, 0,k,l)=A(Lm(ng)-1,Mm(ng)  ,k,l)
                A( 0,-2,k,l)=A(Lm(ng)  ,Mm(ng)-2,k,l)
                A( 0,-1,k,l)=A(Lm(ng)  ,Mm(ng)-1,k,l)
                A( 0, 0,k,l)=A(Lm(ng)  ,Mm(ng)  ,k,l)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE exchange_w4d_tile
#endif

      END MODULE exchange_4d_mod
