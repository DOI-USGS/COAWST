#include "cppdefs.h"
      MODULE ad_exchange_3d_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_exchange_3d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  3D adjoint fields.                                                  !
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
!     ad_A                    3D adjoint field.                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A                    Processed 3D adjoint field.              !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     ad_exchange_p3d_tile    periodic conditions at PSI-points        !
!     ad_exchange_r3d_tile    periodic conditions at RHO-points        !
!     ad_exchange_u3d_tile    periodic conditions at U-points          !
!     ad_exchange_v3d_tile    periodic conditions at V-points          !
!     ad_exchange_w3d_tile    periodic conditions at W-points          !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_exchange_p3d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
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
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,-2,k)=tl_A(Lm(ng)-2,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-2,k)=ad_A(Lm(ng)-2,Mm(ng)-2,k)+      &
     &                                  ad_A(-2,-2,k)
              ad_A(-2,-2,k)=0.0_r8
!>            tl_A(-2,-1,k)=tl_A(Lm(ng)-2,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-1,k)=ad_A(Lm(ng)-2,Mm(ng)-1,k)+      &
     &                                  ad_A(-2,-1,k)
              ad_A(-2,-1,k)=0.0_r8
!>            tl_A(-2, 0,k)=tl_A(Lm(ng)-2,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)  ,k)=ad_A(Lm(ng)-2,Mm(ng)  ,k)+      &
     &                                  ad_A(-2, 0,k)
              ad_A(-2, 0,k)=0.0_r8
!>            tl_A(-1,-2,k)=tl_A(Lm(ng)-1,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-2,k)=ad_A(Lm(ng)-1,Mm(ng)-2,k)+      &
     &                                  ad_A(-1,-2,k)
              ad_A(-1,-2,k)=0.0_r8
!>            tl_A(-1,-1,k)=tl_A(Lm(ng)-1,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-1,k)=ad_A(Lm(ng)-1,Mm(ng)-1,k)+      &
     &                                  ad_A(-1,-1,k)
              ad_A(-1,-1,k)=0.0_r8
!>            tl_A(-1, 0,k)=tl_A(Lm(ng)-1,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)  ,k)=ad_A(Lm(ng)-1,Mm(ng)  ,k)+      &
     &                                  ad_A(-1, 0,k)
              ad_A(-1, 0,k)=0.0_r8
!>            tl_A( 0,-2,k)=tl_A(Lm(ng)  ,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-2,k)=ad_A(Lm(ng)  ,Mm(ng)-2,k)+      &
     &                                  ad_A( 0,-2,k)
              ad_A( 0,-2,k)=0.0_r8
!>            tl_A( 0,-1,k)=tl_A(Lm(ng)  ,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-1,k)=ad_A(Lm(ng)  ,Mm(ng)-1,k)+      &
     &                                  ad_A( 0,-1,k)
              ad_A( 0,-1,k)=0.0_r8
!>            tl_A( 0, 0,k)=tl_A(Lm(ng)  ,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)  ,k)=ad_A(Lm(ng)  ,Mm(ng)  ,k)+      &
     &                                  ad_A( 0, 0,k)
              ad_A( 0, 0,k)=0.0_r8
            END DO
          END IF

          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,-2,k)=tl_A( 1,Mm(ng)-2,k)
!>
              ad_A( 1,Mm(ng)-2,k)=ad_A( 1,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+1,-2,k)
              ad_A(Lm(ng)+1,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+1,-1,k)=tl_A( 1,Mm(ng)-1,k)
!>
              ad_A( 1,Mm(ng)-1,k)=ad_A( 1,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+1,-1,k)
              ad_A(Lm(ng)+1,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+1, 0,k)=tl_A( 1,Mm(ng)  ,k)
!>
              ad_A( 1,Mm(ng)  ,k)=ad_A( 1,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+1, 0,k)
              ad_A(Lm(ng)+1, 0,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-2,k)=tl_A( 2,Mm(ng)-2,k)
!>
              ad_A( 2,Mm(ng)-2,k)=ad_A( 2,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+2,-2,k)
              ad_A(Lm(ng)+2,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-1,k)=tl_A( 2,Mm(ng)-1,k)
!>
              ad_A( 2,Mm(ng)-1,k)=ad_A( 2,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+2,-1,k)
              ad_A(Lm(ng)+2,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+2, 0,k)=tl_A( 2,Mm(ng)  ,k)
!>
              ad_A( 2,Mm(ng)  ,k)=ad_A( 2,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+2, 0,k)
              ad_A(Lm(ng)+2, 0,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+3,-2,k)=tl_A(3 ,Mm(ng)-2,k)
!>
                ad_A(3 ,Mm(ng)-2,k)=ad_A(3 ,Mm(ng)-2,k)+                &
     &                              ad_A(Lm(ng)+3,-2,k)
                ad_A(Lm(ng)+3,-2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,-1,k)=tl_A(3 ,Mm(ng)-1,k)
!>
                ad_A(3 ,Mm(ng)-1,k)=ad_A(3 ,Mm(ng)-1,k)+                &
     &                              ad_A(Lm(ng)+3,-1,k)
                ad_A(Lm(ng)+3,-1,k)=0.0_r8
!>              tl_A(Lm(ng)+3, 0,k)=tl_A(3 ,Mm(ng)  ,k)
!>
                ad_A(3 ,Mm(ng)  ,k)=ad_A(3 ,Mm(ng)  ,k)+                &
     &                              ad_A(Lm(ng)+3, 0,k)
                ad_A(Lm(ng)+3, 0,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,Mm(ng)+1,k)=tl_A(Lm(ng)-2, 1,k)
!>
              ad_A(Lm(ng)-2, 1,k)=ad_A(Lm(ng)-2, 1,k)+                  &
     &                            ad_A(-2,Mm(ng)+1,k)
              ad_A(-2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+1,k)=tl_A(Lm(ng)-1, 1,k)
!>
              ad_A(Lm(ng)-1, 1,k)=ad_A(Lm(ng)-1, 1,k)+                  &
     &                            ad_A(-1,Mm(ng)+1,k)
              ad_A(-1,Mm(ng)+1,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+1,k)=tl_A(Lm(ng)  , 1,k)
!>
              ad_A(Lm(ng)  , 1,k)=ad_A(Lm(ng)  , 1,k)+                  &
     &                            ad_A( 0,Mm(ng)+1,k)
              ad_A( 0,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-2,Mm(ng)+2,k)=tl_A(Lm(ng)-2, 2,k)
!>
              ad_A(Lm(ng)-2, 2,k)=ad_A(Lm(ng)-2, 2,k)+                  &
     &                            ad_A(-2,Mm(ng)+2,k)
              ad_A(-2,Mm(ng)+2,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+2,k)=tl_A(Lm(ng)-1, 2,k)
!>
              ad_A(Lm(ng)-1, 2,k)=ad_A(Lm(ng)-1, 2,k)+                  &
     &                            ad_A(-1,Mm(ng)+2,k)
              ad_A(-1,Mm(ng)+2,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+2,k)=tl_A(Lm(ng)  , 2,k)
!>
              ad_A(Lm(ng)  , 2,k)=ad_A(Lm(ng)  , 2,k)+                  &
     &                            ad_A( 0,Mm(ng)+2,k)
              ad_A( 0,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(-2,Mm(ng)+3,k)=tl_A(Lm(ng)-2, 3,k)
!>
                ad_A(Lm(ng)-2, 3,k)=ad_A(Lm(ng)-2, 3,k)+                &
     &                              ad_A(-2,Mm(ng)+3,k)
                ad_A(-2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(-1,Mm(ng)+3,k)=tl_A(Lm(ng)-1, 3,k)
!>
                ad_A(Lm(ng)-1, 3,k)=ad_A(Lm(ng)-1, 3,k)+                &
     &                              ad_A(-1,Mm(ng)+3,k)
                ad_A(-1,Mm(ng)+3,k)=0.0_r8
!>              tl_A( 0,Mm(ng)+3,k)=tl_A(Lm(ng)  , 3,k)
!>
                ad_A(Lm(ng)  , 3,k)=ad_A(Lm(ng)  , 3,k)+                &
     &                              ad_A( 0,Mm(ng)+3,k)
                ad_A( 0,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,Mm(ng)+1,k)=tl_A( 1, 1,k)
!>
              ad_A( 1, 1,k)=ad_A( 1, 1,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+1,k)
              ad_A(Lm(ng)+1,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+1,Mm(ng)+2,k)=tl_A( 1, 2,k)
!>
              ad_A( 1, 2,k)=ad_A( 1, 2,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+2,k)
              ad_A(Lm(ng)+1,Mm(ng)+2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+1,k)=tl_A( 2, 1,k)
!>
              ad_A( 2, 1,k)=ad_A( 2, 1,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+1,k)
              ad_A(Lm(ng)+2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+2,k)=tl_A(2, 2,k)
!>
              ad_A( 2, 2,k)=ad_A( 2, 2,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+2,k)
              ad_A(Lm(ng)+2,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+1,Mm(ng)+3,k)=tl_A( 1, 3,k)
!>
                ad_A( 1, 3,k)=ad_A( 1, 3,k)+                            &
     &                        ad_A(Lm(ng)+1,Mm(ng)+3,k)
                ad_A(Lm(ng)+1,Mm(ng)+3,k)=0.0_r8

!>              tl_A(Lm(ng)+2,Mm(ng)+3,k)=tl_A( 2, 3,k)
!>
                ad_A( 2, 3,k)=ad_A( 2, 3,k)+                            &
     &                        ad_A(Lm(ng)+2,Mm(ng)+3,k)
                ad_A(Lm(ng)+2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+1,k)=tl_A( 3, 1,k)
!>
                ad_A( 3, 1,k)=ad_A( 3, 1,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+1,k)
                ad_A(Lm(ng)+3,Mm(ng)+1,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+2,k)=tl_A( 3, 2,k)
!>
                ad_A( 3, 2,k)=ad_A( 3, 2,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+2,k)
                ad_A(Lm(ng)+3,Mm(ng)+2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+3,k)=tl_A( 3, 3,k)
!>
                ad_A( 3, 3,k)=ad_A( 3, 3,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+3,k)
                ad_A(Lm(ng)+3,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
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
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,-2,k)=tl_A(i,Mm(ng)-2,k)
!>
                ad_A(i,Mm(ng)-2,k)=ad_A(i,Mm(ng)-2,k)+                  &
     &                             ad_A(i,-2,k)
                ad_A(i,-2,k)=0.0_r8
!>              tl_A(i,-1,k)=tl_A(i,Mm(ng)-1,k)
!>
                ad_A(i,Mm(ng)-1,k)=ad_A(i,Mm(ng)-1,k)+                  &
     &                             ad_A(i,-1,k)
                ad_A(i,-1,k)=0.0_r8
!>              tl_A(i, 0,k)=tl_A(i,Mm(ng)  ,k)
!>
                ad_A(i,Mm(ng)  ,k)=ad_A(i,Mm(ng)  ,k)+                  &
     &                             ad_A(i, 0,k)
                ad_A(i, 0,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,Mm(ng)+1,k)=tl_A(i,1,k)
!>
                ad_A(i,1,k)=ad_A(i,1,k)+                                &
     &                      ad_A(i,Mm(ng)+1,k)
                ad_A(i,Mm(ng)+1,k)=0.0_r8
!>              tl_A(i,Mm(ng)+2,k)=tl_A(i,2,k)
!>
                ad_A(i,2,k)=ad_A(i,2,k)+                                &
     &                      ad_A(i,Mm(ng)+2,k)
                ad_A(i,Mm(ng)+2,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax

!>                tl_A(i,Mm(ng)+3,k)=tl_A(i,3,k)
!>
                  ad_A(i,3,k)=ad_A(i,3,k)+                              &
     &                        ad_A(i,Mm(ng)+3,k)
                  ad_A(i,Mm(ng)+3,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
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
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(-2,j,k)=tl_A(Lm(ng)-2,j,k)
!>
                ad_A(Lm(ng)-2,j,k)=ad_A(Lm(ng)-2,j,k)+                  &
     &                             ad_A(-2,j,k)
                ad_A(-2,j,k)=0.0_r8
!>              tl_A(-1,j,k)=tl_A(Lm(ng)-1,j,k)
!>
                ad_A(Lm(ng)-1,j,k)=ad_A(Lm(ng)-1,j,k)+                  &
     &                             ad_A(-1,j,k)
                ad_A(-1,j,k)=0.0_r8
!>              tl_A( 0,j,k)=tl_A(Lm(ng)  ,j,k)
!>
                ad_A(Lm(ng)  ,j,k)=ad_A(Lm(ng)  ,j,k)+                  &
     &                             ad_A( 0,j,k)
                ad_A( 0,j,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(Lm(ng)+1,j,k)=tl_A(1,j,k)
!>
                ad_A(1,j,k)=ad_A(1,j,k)+                                &
     &                      ad_A(Lm(ng)+1,j,k)
                ad_A(Lm(ng)+1,j,k)=0.0_r8
!>              tl_A(Lm(ng)+2,j,k)=tl_A(2,j,k)
!>
                ad_A(2,j,k)=ad_A(2,j,k)+                                &
     &                      ad_A(Lm(ng)+2,j,k)
                ad_A(Lm(ng)+2,j,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
!>                tl_A(Lm(ng)+3,j,k)=tl_A(3,j,k)
!>
                  ad_A(3,j,k)=ad_A(3,j,k)+                              &
     &                        ad_A(Lm(ng)+3,j,k)
                  ad_A(Lm(ng)+3,j,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_exchange_p3d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_r3d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
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
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,-2,k)=tl_A(Lm(ng)-2,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-2,k)=ad_A(Lm(ng)-2,Mm(ng)-2,k)+      &
     &                                  ad_A(-2,-2,k)
              ad_A(-2,-2,k)=0.0_r8
!>            tl_A(-2,-1,k)=tl_A(Lm(ng)-2,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-1,k)=ad_A(Lm(ng)-2,Mm(ng)-1,k)+      &
     &                                  ad_A(-2,-1,k)
              ad_A(-2,-1,k)=0.0_r8
!>            tl_A(-2, 0,k)=tl_A(Lm(ng)-2,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)  ,k)=ad_A(Lm(ng)-2,Mm(ng)  ,k)+      &
     &                                  ad_A(-2, 0,k)
              ad_A(-2, 0,k)=0.0_r8
!>            tl_A(-1,-2,k)=tl_A(Lm(ng)-1,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-2,k)=ad_A(Lm(ng)-1,Mm(ng)-2,k)+      &
     &                                  ad_A(-1,-2,k)
              ad_A(-1,-2,k)=0.0_r8
!>            tl_A(-1,-1,k)=tl_A(Lm(ng)-1,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-1,k)=ad_A(Lm(ng)-1,Mm(ng)-1,k)+      &
     &                                  ad_A(-1,-1,k)
              ad_A(-1,-1,k)=0.0_r8
!>            tl_A(-1, 0,k)=tl_A(Lm(ng)-1,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)  ,k)=ad_A(Lm(ng)-1,Mm(ng)  ,k)+      &
     &                                  ad_A(-1, 0,k)
              ad_A(-1, 0,k)=0.0_r8
!>            tl_A( 0,-2,k)=tl_A(Lm(ng)  ,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-2,k)=ad_A(Lm(ng)  ,Mm(ng)-2,k)+      &
     &                                  ad_A( 0,-2,k)
              ad_A( 0,-2,k)=0.0_r8
!>            tl_A( 0,-1,k)=tl_A(Lm(ng)  ,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-1,k)=ad_A(Lm(ng)  ,Mm(ng)-1,k)+      &
     &                                  ad_A( 0,-1,k)
              ad_A( 0,-1,k)=0.0_r8
!>            tl_A( 0, 0,k)=tl_A(Lm(ng)  ,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)  ,k)=ad_A(Lm(ng)  ,Mm(ng)  ,k)+      &
     &                                  ad_A( 0, 0,k)
              ad_A( 0, 0,k)=0.0_r8
            END DO
          END IF

          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,-2,k)=tl_A( 1,Mm(ng)-2,k)
!>
              ad_A( 1,Mm(ng)-2,k)=ad_A( 1,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+1,-2,k)
              ad_A(Lm(ng)+1,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+1,-1,k)=tl_A( 1,Mm(ng)-1,k)
!>
              ad_A( 1,Mm(ng)-1,k)=ad_A( 1,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+1,-1,k)
              ad_A(Lm(ng)+1,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+1, 0,k)=tl_A( 1,Mm(ng)  ,k)
!>
              ad_A( 1,Mm(ng)  ,k)=ad_A( 1,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+1, 0,k)
              ad_A(Lm(ng)+1, 0,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-2,k)=tl_A( 2,Mm(ng)-2,k)
!>
              ad_A( 2,Mm(ng)-2,k)=ad_A( 2,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+2,-2,k)
              ad_A(Lm(ng)+2,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-1,k)=tl_A( 2,Mm(ng)-1,k)
!>
              ad_A( 2,Mm(ng)-1,k)=ad_A( 2,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+2,-1,k)
              ad_A(Lm(ng)+2,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+2, 0,k)=tl_A( 2,Mm(ng)  ,k)
!>
              ad_A( 2,Mm(ng)  ,k)=ad_A( 2,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+2, 0,k)
              ad_A(Lm(ng)+2, 0,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+3,-2,k)=tl_A(3 ,Mm(ng)-2,k)
!>
                ad_A(3 ,Mm(ng)-2,k)=ad_A(3 ,Mm(ng)-2,k)+                &
     &                              ad_A(Lm(ng)+3,-2,k)
                ad_A(Lm(ng)+3,-2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,-1,k)=tl_A(3 ,Mm(ng)-1,k)
!>
                ad_A(3 ,Mm(ng)-1,k)=ad_A(3 ,Mm(ng)-1,k)+                &
     &                              ad_A(Lm(ng)+3,-1,k)
                ad_A(Lm(ng)+3,-1,k)=0.0_r8
!>              tl_A(Lm(ng)+3, 0,k)=tl_A(3 ,Mm(ng)  ,k)
!>
                ad_A(3 ,Mm(ng)  ,k)=ad_A(3 ,Mm(ng)  ,k)+                &
     &                              ad_A(Lm(ng)+3, 0,k)
                ad_A(Lm(ng)+3, 0,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,Mm(ng)+1,k)=tl_A(Lm(ng)-2, 1,k)
!>
              ad_A(Lm(ng)-2, 1,k)=ad_A(Lm(ng)-2, 1,k)+                  &
     &                            ad_A(-2,Mm(ng)+1,k)
              ad_A(-2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+1,k)=tl_A(Lm(ng)-1, 1,k)
!>
              ad_A(Lm(ng)-1, 1,k)=ad_A(Lm(ng)-1, 1,k)+                  &
     &                            ad_A(-1,Mm(ng)+1,k)
              ad_A(-1,Mm(ng)+1,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+1,k)=tl_A(Lm(ng)  , 1,k)
!>
              ad_A(Lm(ng)  , 1,k)=ad_A(Lm(ng)  , 1,k)+                  &
     &                            ad_A( 0,Mm(ng)+1,k)
              ad_A( 0,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-2,Mm(ng)+2,k)=tl_A(Lm(ng)-2, 2,k)
!>
              ad_A(Lm(ng)-2, 2,k)=ad_A(Lm(ng)-2, 2,k)+                  &
     &                            ad_A(-2,Mm(ng)+2,k)
              ad_A(-2,Mm(ng)+2,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+2,k)=tl_A(Lm(ng)-1, 2,k)
!>
              ad_A(Lm(ng)-1, 2,k)=ad_A(Lm(ng)-1, 2,k)+                  &
     &                            ad_A(-1,Mm(ng)+2,k)
              ad_A(-1,Mm(ng)+2,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+2,k)=tl_A(Lm(ng)  , 2,k)
!>
              ad_A(Lm(ng)  , 2,k)=ad_A(Lm(ng)  , 2,k)+                  &
     &                            ad_A( 0,Mm(ng)+2,k)
              ad_A( 0,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(-2,Mm(ng)+3,k)=tl_A(Lm(ng)-2, 3,k)
!>
                ad_A(Lm(ng)-2, 3,k)=ad_A(Lm(ng)-2, 3,k)+                &
     &                              ad_A(-2,Mm(ng)+3,k)
                ad_A(-2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(-1,Mm(ng)+3,k)=tl_A(Lm(ng)-1, 3,k)
!>
                ad_A(Lm(ng)-1, 3,k)=ad_A(Lm(ng)-1, 3,k)+                &
     &                              ad_A(-1,Mm(ng)+3,k)
                ad_A(-1,Mm(ng)+3,k)=0.0_r8
!>              tl_A( 0,Mm(ng)+3,k)=tl_A(Lm(ng)  , 3,k)
!>
                ad_A(Lm(ng)  , 3,k)=ad_A(Lm(ng)  , 3,k)+                &
     &                              ad_A( 0,Mm(ng)+3,k)
                ad_A( 0,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,Mm(ng)+1,k)=tl_A( 1, 1,k)
!>
              ad_A( 1, 1,k)=ad_A( 1, 1,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+1,k)
              ad_A(Lm(ng)+1,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+1,Mm(ng)+2,k)=tl_A( 1, 2,k)
!>
              ad_A( 1, 2,k)=ad_A( 1, 2,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+2,k)
              ad_A(Lm(ng)+1,Mm(ng)+2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+1,k)=tl_A( 2, 1,k)
!>
              ad_A( 2, 1,k)=ad_A( 2, 1,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+1,k)
              ad_A(Lm(ng)+2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+2,k)=tl_A( 2, 2,k)
!>
              ad_A( 2, 2,k)=ad_A( 2, 2,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+2,k)
              ad_A(Lm(ng)+2,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+1,Mm(ng)+3,k)=tl_A( 1, 3,k)
!>
                ad_A( 1, 3,k)=ad_A( 1, 3,k)+                            &
     &                        ad_A(Lm(ng)+1,Mm(ng)+3,k)
                ad_A(Lm(ng)+1,Mm(ng)+3,k)=0.0_r8

!>              tl_A(Lm(ng)+2,Mm(ng)+3,k)=tl_A( 2, 3,k)
!>
                ad_A( 2, 3,k)=ad_A( 2, 3,k)+                            &
     &                        ad_A(Lm(ng)+2,Mm(ng)+3,k)
                ad_A(Lm(ng)+2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+1,k)=tl_A( 3, 1,k)
!>
                ad_A( 3, 1,k)=ad_A( 3, 1,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+1,k)
                ad_A(Lm(ng)+3,Mm(ng)+1,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+2,k)=tl_A( 3, 2,k)
!>
                ad_A( 3, 2,k)=ad_A( 3, 2,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+2,k)
                ad_A(Lm(ng)+3,Mm(ng)+2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+3,k)=tl_A( 3, 3,k)
!>
                ad_A( 3, 3,k)=ad_A( 3, 3,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+3,k)
                ad_A(Lm(ng)+3,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
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
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,-2,k)=tl_A(i,Mm(ng)-2,k)
!>
                ad_A(i,Mm(ng)-2,k)=ad_A(i,Mm(ng)-2,k)+                  &
     &                             ad_A(i,-2,k)
                ad_A(i,-2,k)=0.0_r8
!>              tl_A(i,-1,k)=tl_A(i,Mm(ng)-1,k)
!>
                ad_A(i,Mm(ng)-1,k)=ad_A(i,Mm(ng)-1,k)+                  &
     &                             ad_A(i,-1,k)
                ad_A(i,-1,k)=0.0_r8
!>              tl_A(i, 0,k)=tl_A(i,Mm(ng)  ,k)
!>
                ad_A(i,Mm(ng)  ,k)=ad_A(i,Mm(ng)  ,k)+                  &
     &                             ad_A(i, 0,k)
                ad_A(i, 0,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,Mm(ng)+1,k)=tl_A(i, 1,k)
!>
                ad_A(i, 1,k)=ad_A(i, 1,k)+                              &
     &                       ad_A(i,Mm(ng)+1,k)
                ad_A(i,Mm(ng)+1,k)=0.0_r8
!>              tl_A(i,Mm(ng)+2,k)=tl_A(i, 2,k)
!>
                ad_A(i, 2,k)=ad_A(i, 2,k)+                              &
     &                       ad_A(i,Mm(ng)+2,k)
                ad_A(i,Mm(ng)+2,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
!>                tl_A(i,Mm(ng)+3,k)=tl_A(i,3,k)
!>
                  ad_A(i,3,k)=ad_A(i,3,k)+                              &
     &                        ad_A(i,Mm(ng)+3,k)
                  ad_A(i,Mm(ng)+3,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
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
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(-2,j,k)=tl_A(Lm(ng)-2,j,k)
!>
                ad_A(Lm(ng)-2,j,k)=ad_A(Lm(ng)-2,j,k)+                  &
     &                             ad_A(-2,j,k)
                ad_A(-2,j,k)=0.0_r8
!>              tl_A(-1,j,k)=tl_A(Lm(ng)-1,j,k)
!>
                ad_A(Lm(ng)-1,j,k)=ad_A(Lm(ng)-1,j,k)+                  &
     &                             ad_A(-1,j,k)
                ad_A(-1,j,k)=0.0_r8
!>              tl_A( 0,j,k)=tl_A(Lm(ng)  ,j,k)
!>
                ad_A(Lm(ng)  ,j,k)=ad_A(Lm(ng)  ,j,k)+                  &
     &                             ad_A( 0,j,k)
                ad_A( 0,j,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(Lm(ng)+1,j,k)=tl_A(1,j,k)
!>
                ad_A(1,j,k)=ad_A(1,j,k)+                                &
     &                      ad_A(Lm(ng)+1,j,k)
                ad_A(Lm(ng)+1,j,k)=0.0_r8
!>              tl_A(Lm(ng)+2,j,k)=tl_A( 2,j,k)
!>
                ad_A(2,j,k)=ad_A(2,j,k)+                                &
     &                      ad_A(Lm(ng)+2,j,k)
                ad_A(Lm(ng)+2,j,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
!>                tl_A(Lm(ng)+3,j,k)=tl_A(3,j,k)
!>
                  ad_A(3,j,k)=ad_A(3,j,k)+                              &
     &                        ad_A(Lm(ng)+3,j,k)
                  ad_A(Lm(ng)+3,j,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_exchange_r3d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_u3d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
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
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,-2,k)=tl_A(Lm(ng)-2,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-2,k)=ad_A(Lm(ng)-2,Mm(ng)-2,k)+      &
     &                                  ad_A(-2,-2,k)
              ad_A(-2,-2,k)=0.0_r8
!>            tl_A(-2,-1,k)=tl_A(Lm(ng)-2,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-1,k)=ad_A(Lm(ng)-2,Mm(ng)-1,k)+      &
     &                                  ad_A(-2,-1,k)
              ad_A(-2,-1,k)=0.0_r8
!>            tl_A(-2, 0,k)=tl_A(Lm(ng)-2,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)  ,k)=ad_A(Lm(ng)-2,Mm(ng)  ,k)+      &
     &                                  ad_A(-2, 0,k)
              ad_A(-2, 0,k)=0.0_r8
!>            tl_A(-1,-2,k)=tl_A(Lm(ng)-1,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-2,k)=ad_A(Lm(ng)-1,Mm(ng)-2,k)+      &
     &                                  ad_A(-1,-2,k)
              ad_A(-1,-2,k)=0.0_r8
!>            tl_A(-1,-1,k)=tl_A(Lm(ng)-1,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-1,k)=ad_A(Lm(ng)-1,Mm(ng)-1,k)+      &
     &                                  ad_A(-1,-1,k)
              ad_A(-1,-1,k)=0.0_r8
!>            tl_A(-1, 0,k)=tl_A(Lm(ng)-1,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)  ,k)=ad_A(Lm(ng)-1,Mm(ng)  ,k)+      &
     &                                  ad_A(-1, 0,k)
              ad_A(-1, 0,k)=0.0_r8
!>            tl_A( 0,-2,k)=tl_A(Lm(ng)  ,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-2,k)=ad_A(Lm(ng)  ,Mm(ng)-2,k)+      &
     &                                  ad_A( 0,-2,k)
              ad_A( 0,-2,k)=0.0_r8
!>            tl_A( 0,-1,k)=tl_A(Lm(ng)  ,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-1,k)=ad_A(Lm(ng)  ,Mm(ng)-1,k)+      &
     &                                  ad_A( 0,-1,k)
              ad_A( 0,-1,k)=0.0_r8
!>            tl_A( 0, 0,k)=tl_A(Lm(ng)  ,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)  ,k)=ad_A(Lm(ng)  ,Mm(ng)  ,k)+      &
     &                                  ad_A( 0, 0,k)
              ad_A( 0, 0,k)=0.0_r8
            END DO
          END IF

          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,-2,k)=tl_A( 1,Mm(ng)-2,k)
!>
              ad_A( 1,Mm(ng)-2,k)=ad_A( 1,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+1,-2,k)
              ad_A(Lm(ng)+1,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+1,-1,k)=tl_A( 1,Mm(ng)-1,k)
!>
              ad_A( 1,Mm(ng)-1,k)=ad_A( 1,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+1,-1,k)
              ad_A(Lm(ng)+1,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+1, 0,k)=tl_A( 1,Mm(ng)  ,k)
!>
              ad_A( 1,Mm(ng)  ,k)=ad_A( 1,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+1, 0,k)
              ad_A(Lm(ng)+1, 0,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-2,k)=tl_A( 2,Mm(ng)-2,k)
!>
              ad_A( 2,Mm(ng)-2,k)=ad_A( 2,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+2,-2,k)
              ad_A(Lm(ng)+2,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-1,k)=tl_A( 2,Mm(ng)-1,k)
!>
              ad_A( 2,Mm(ng)-1,k)=ad_A( 2,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+2,-1,k)
              ad_A(Lm(ng)+2,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+2, 0,k)=tl_A( 2,Mm(ng)  ,k)
!>
              ad_A( 2,Mm(ng)  ,k)=ad_A( 2,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+2, 0,k)
              ad_A(Lm(ng)+2, 0,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+3,-2,k)=tl_A(3 ,Mm(ng)-2,k)
!>
                ad_A(3 ,Mm(ng)-2,k)=ad_A(3 ,Mm(ng)-2,k)+                &
     &                              ad_A(Lm(ng)+3,-2,k)
                ad_A(Lm(ng)+3,-2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,-1,k)=tl_A(3 ,Mm(ng)-1,k)
!>
                ad_A(3 ,Mm(ng)-1,k)=ad_A(3 ,Mm(ng)-1,k)+                &
     &                              ad_A(Lm(ng)+3,-1,k)
                ad_A(Lm(ng)+3,-1,k)=0.0_r8
!>              tl_A(Lm(ng)+3, 0,k)=tl_A(3 ,Mm(ng)  ,k)
!>
                ad_A(3 ,Mm(ng)  ,k)=ad_A(3 ,Mm(ng)  ,k)+                &
     &                              ad_A(Lm(ng)+3, 0,k)
                ad_A(Lm(ng)+3, 0,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,Mm(ng)+1,k)=tl_A(Lm(ng)-2, 1,k)
!>
              ad_A(Lm(ng)-2, 1,k)=ad_A(Lm(ng)-2, 1,k)+                  &
     &                            ad_A(-2,Mm(ng)+1,k)
              ad_A(-2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+1,k)=tl_A(Lm(ng)-1, 1,k)
!>
              ad_A(Lm(ng)-1, 1,k)=ad_A(Lm(ng)-1, 1,k)+                  &
     &                            ad_A(-1,Mm(ng)+1,k)
              ad_A(-1,Mm(ng)+1,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+1,k)=tl_A(Lm(ng)  , 1,k)
!>
              ad_A(Lm(ng)  , 1,k)=ad_A(Lm(ng)  , 1,k)+                  &
     &                            ad_A( 0,Mm(ng)+1,k)
              ad_A( 0,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-2,Mm(ng)+2,k)=tl_A(Lm(ng)-2, 2,k)
!>
              ad_A(Lm(ng)-2, 2,k)=ad_A(Lm(ng)-2, 2,k)+                  &
     &                            ad_A(-2,Mm(ng)+2,k)
              ad_A(-2,Mm(ng)+2,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+2,k)=tl_A(Lm(ng)-1, 2,k)
!>
              ad_A(Lm(ng)-1, 2,k)=ad_A(Lm(ng)-1, 2,k)+                  &
     &                            ad_A(-1,Mm(ng)+2,k)
              ad_A(-1,Mm(ng)+2,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+2,k)=tl_A(Lm(ng)  , 2,k)
!>
              ad_A(Lm(ng)  , 2,k)=ad_A(Lm(ng)  , 2,k)+                  &
     &                            ad_A( 0,Mm(ng)+2,k)
              ad_A( 0,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(-2,Mm(ng)+3,k)=tl_A(Lm(ng)-2, 3,k)
!>
                ad_A(Lm(ng)-2, 3,k)=ad_A(Lm(ng)-2, 3,k)+                &
     &                              ad_A(-2,Mm(ng)+3,k)
                ad_A(-2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(-1,Mm(ng)+3,k)=tl_A(Lm(ng)-1, 3,k)
!>
                ad_A(Lm(ng)-1, 3,k)=ad_A(Lm(ng)-1, 3,k)+                &
     &                              ad_A(-1,Mm(ng)+3,k)
                ad_A(-1,Mm(ng)+3,k)=0.0_r8
!>              tl_A( 0,Mm(ng)+3,k)=tl_A(Lm(ng)  , 3,k)
!>
                ad_A(Lm(ng)  , 3,k)=ad_A(Lm(ng)  , 3,k)+                &
     &                              ad_A( 0,Mm(ng)+3,k)
                ad_A( 0,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,Mm(ng)+1,k)=tl_A( 1  , 1  ,k)
!>
              ad_A( 1, 1,k)=ad_A( 1, 1,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+1,k)
              ad_A(Lm(ng)+1,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+1,Mm(ng)+2,k)=tl_A( 1  , 2,k)
!>
              ad_A( 1, 2,k)=ad_A( 1, 2,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+2,k)
              ad_A(Lm(ng)+1,Mm(ng)+2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+1,k)=tl_A( 2, 1  ,k)
!>
              ad_A( 2, 1,k)=ad_A( 2, 1,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+1,k)
              ad_A(Lm(ng)+2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+2,k)=tl_A( 2, 2,k)
!>
              ad_A( 2, 2,k)=ad_A( 2, 2,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+2,k)
              ad_A(Lm(ng)+2,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+1,Mm(ng)+3,k)=tl_A( 1, 3,k)
!>
                ad_A( 1, 3,k)=ad_A( 1, 3,k)+                            &
     &                        ad_A(Lm(ng)+1,Mm(ng)+3,k)
                ad_A(Lm(ng)+1,Mm(ng)+3,k)=0.0_r8

!>              tl_A(Lm(ng)+2,Mm(ng)+3,k)=tl_A( 2, 3,k)
!>
                ad_A( 2, 3,k)=ad_A( 2, 3,k)+                            &
     &                        ad_A(Lm(ng)+2,Mm(ng)+3,k)
                ad_A(Lm(ng)+2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+1,k)=tl_A( 3, 1,k)
!>
                ad_A( 3, 1,k)=ad_A( 3, 1,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+1,k)
                ad_A(Lm(ng)+3,Mm(ng)+1,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+2,k)=tl_A( 3, 2,k)
!>
                ad_A( 3, 2,k)=ad_A( 3, 2,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+2,k)
                ad_A(Lm(ng)+3,Mm(ng)+2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+3,k)=tl_A( 3, 3,k)
!>
                ad_A( 3, 3,k)=ad_A( 3, 3,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+3,k)
                ad_A(Lm(ng)+3,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
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
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,-2,k)=tl_A(i,Mm(ng)-2,k)
!>
                ad_A(i,Mm(ng)-2,k)=ad_A(i,Mm(ng)-2,k)+                  &
     &                             ad_A(i,-2,k)
                ad_A(i,-2,k)=0.0_r8
!>              tl_A(i,-1,k)=tl_A(i,Mm(ng)-1,k)
!>
                ad_A(i,Mm(ng)-1,k)=ad_A(i,Mm(ng)-1,k)+                  &
     &                             ad_A(i,-1,k)
                ad_A(i,-1,k)=0.0_r8
!>              tl_A(i, 0,k)=tl_A(i,Mm(ng)  ,k)
!>
                ad_A(i,Mm(ng)  ,k)=ad_A(i,Mm(ng)  ,k)+                  &
     &                             ad_A(i, 0,k)
                ad_A(i, 0,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,Mm(ng)+1,k)=tl_A(i,1,k)
!>
                ad_A(i,1,k)=ad_A(i,1,k)+                                &
     &                      ad_A(i,Mm(ng)+1,k)
                ad_A(i,Mm(ng)+1,k)=0.0_r8
!>              tl_A(i,Mm(ng)+2,k)=tl_A(i, 2,k)
!>
                ad_A(i,2,k)=ad_A(i,2,k)+                                &
     &                      ad_A(i,Mm(ng)+2,k)
                ad_A(i,Mm(ng)+2,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
!>                tl_A(i,Mm(ng)+3,k)=tl_A(i,3,k)
!>
                  ad_A(i,3,k)=ad_A(i,3,k)+                              &
     &                        ad_A(i,Mm(ng)+3,k)
                  ad_A(i,Mm(ng)+3,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
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
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(-2,j,k)=tl_A(Lm(ng)-2,j,k)
!>
                ad_A(Lm(ng)-2,j,k)=ad_A(Lm(ng)-2,j,k)+                  &
     &                             ad_A(-2,j,k)
                ad_A(-2,j,k)=0.0_r8
!>              tl_A(-1,j,k)=tl_A(Lm(ng)-1,j,k)
!>
                ad_A(Lm(ng)-1,j,k)=ad_A(Lm(ng)-1,j,k)+                  &
     &                             ad_A(-1,j,k)
                ad_A(-1,j,k)=0.0_r8
!>              tl_A( 0,j,k)=tl_A(Lm(ng)  ,j,k)
!>
                ad_A(Lm(ng)  ,j,k)=ad_A(Lm(ng)  ,j,k)+                  &
     &                             ad_A( 0,j,k)
                ad_A( 0,j,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(Lm(ng)+1,j,k)=tl_A(1,j,k)
!>
                ad_A(1,j,k)=ad_A(1,j,k)+                                &
     &                      ad_A(Lm(ng)+1,j,k)
                ad_A(Lm(ng)+1,j,k)=0.0_r8
!>              tl_A(Lm(ng)+2,j,k)=tl_A( 2,j,k)
!>
                ad_A(2,j,k)=ad_A(2,j,k)+                                &
     &                      ad_A(Lm(ng)+2,j,k)
                ad_A(Lm(ng)+2,j,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
!>                tl_A(Lm(ng)+3,j,k)=tl_A(3,j,k)
!>
                  ad_A(3,j,k)=ad_A(3,j,k)+                              &
     &                        ad_A(Lm(ng)+3,j,k)
                  ad_A(Lm(ng)+3,j,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_exchange_u3d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_v3d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
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
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,-2,k)=tl_A(Lm(ng)-2,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-2,k)=ad_A(Lm(ng)-2,Mm(ng)-2,k)+      &
     &                                  ad_A(-2,-2,k)
              ad_A(-2,-2,k)=0.0_r8
!>            tl_A(-2,-1,k)=tl_A(Lm(ng)-2,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-1,k)=ad_A(Lm(ng)-2,Mm(ng)-1,k)+      &
     &                                  ad_A(-2,-1,k)
              ad_A(-2,-1,k)=0.0_r8
!>            tl_A(-2, 0,k)=tl_A(Lm(ng)-2,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)  ,k)=ad_A(Lm(ng)-2,Mm(ng)  ,k)+      &
     &                                  ad_A(-2, 0,k)
              ad_A(-2, 0,k)=0.0_r8
!>            tl_A(-1,-2,k)=tl_A(Lm(ng)-1,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-2,k)=ad_A(Lm(ng)-1,Mm(ng)-2,k)+      &
     &                                  ad_A(-1,-2,k)
              ad_A(-1,-2,k)=0.0_r8
!>            tl_A(-1,-1,k)=tl_A(Lm(ng)-1,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-1,k)=ad_A(Lm(ng)-1,Mm(ng)-1,k)+      &
     &                                  ad_A(-1,-1,k)
              ad_A(-1,-1,k)=0.0_r8
!>            tl_A(-1, 0,k)=tl_A(Lm(ng)-1,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)  ,k)=ad_A(Lm(ng)-1,Mm(ng)  ,k)+      &
     &                                  ad_A(-1, 0,k)
              ad_A(-1, 0,k)=0.0_r8
!>            tl_A( 0,-2,k)=tl_A(Lm(ng)  ,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-2,k)=ad_A(Lm(ng)  ,Mm(ng)-2,k)+      &
     &                                  ad_A( 0,-2,k)
              ad_A( 0,-2,k)=0.0_r8
!>            tl_A( 0,-1,k)=tl_A(Lm(ng)  ,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-1,k)=ad_A(Lm(ng)  ,Mm(ng)-1,k)+      &
     &                                  ad_A( 0,-1,k)
              ad_A( 0,-1,k)=0.0_r8
!>            tl_A( 0, 0,k)=tl_A(Lm(ng)  ,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)  ,k)=ad_A(Lm(ng)  ,Mm(ng)  ,k)+      &
     &                                  ad_A( 0, 0,k)
              ad_A( 0, 0,k)=0.0_r8
            END DO
          END IF

          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,-2,k)=tl_A( 1,Mm(ng)-2,k)
!>
              ad_A( 1,Mm(ng)-2,k)=ad_A( 1,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+1,-2,k)
              ad_A(Lm(ng)+1,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+1,-1,k)=tl_A( 1,Mm(ng)-1,k)
!>
              ad_A( 1,Mm(ng)-1,k)=ad_A( 1,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+1,-1,k)
              ad_A(Lm(ng)+1,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+1, 0,k)=tl_A( 1,Mm(ng)  ,k)
!>
              ad_A( 1,Mm(ng)  ,k)=ad_A( 1,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+1, 0,k)
              ad_A(Lm(ng)+1, 0,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-2,k)=tl_A( 2,Mm(ng)-2,k)
!>
              ad_A( 2,Mm(ng)-2,k)=ad_A( 2,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+2,-2,k)
              ad_A(Lm(ng)+2,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-1,k)=tl_A( 2,Mm(ng)-1,k)
!>
              ad_A( 2,Mm(ng)-1,k)=ad_A( 2,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+2,-1,k)
              ad_A(Lm(ng)+2,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+2, 0,k)=tl_A( 2,Mm(ng)  ,k)
!>
              ad_A( 2,Mm(ng)  ,k)=ad_A( 2,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+2, 0,k)
              ad_A(Lm(ng)+2, 0,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+3,-2,k)=tl_A(3 ,Mm(ng)-2,k)
!>
                ad_A(3 ,Mm(ng)-2,k)=ad_A(3 ,Mm(ng)-2,k)+                &
     &                              ad_A(Lm(ng)+3,-2,k)
                ad_A(Lm(ng)+3,-2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,-1,k)=tl_A(3 ,Mm(ng)-1,k)
!>
                ad_A(3 ,Mm(ng)-1,k)=ad_A(3 ,Mm(ng)-1,k)+                &
     &                              ad_A(Lm(ng)+3,-1,k)
                ad_A(Lm(ng)+3,-1,k)=0.0_r8
!>              tl_A(Lm(ng)+3, 0,k)=tl_A(3 ,Mm(ng)  ,k)
!>
                ad_A(3 ,Mm(ng)  ,k)=ad_A(3 ,Mm(ng)  ,k)+                &
     &                              ad_A(Lm(ng)+3, 0,k)
                ad_A(Lm(ng)+3, 0,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,Mm(ng)+1,k)=tl_A(Lm(ng)-2, 1,k)
!>
              ad_A(Lm(ng)-2, 1,k)=ad_A(Lm(ng)-2, 1,k)+                  &
     &                            ad_A(-2,Mm(ng)+1,k)
              ad_A(-2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+1,k)=tl_A(Lm(ng)-1, 1,k)
!>
              ad_A(Lm(ng)-1, 1,k)=ad_A(Lm(ng)-1, 1,k)+                  &
     &                            ad_A(-1,Mm(ng)+1,k)
              ad_A(-1,Mm(ng)+1,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+1,k)=tl_A(Lm(ng)  , 1,k)
!>
              ad_A(Lm(ng)  , 1,k)=ad_A(Lm(ng)  , 1,k)+                  &
     &                            ad_A( 0,Mm(ng)+1,k)
              ad_A( 0,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-2,Mm(ng)+2,k)=tl_A(Lm(ng)-2, 2,k)
!>
              ad_A(Lm(ng)-2, 2,k)=ad_A(Lm(ng)-2, 2,k)+                  &
     &                            ad_A(-2,Mm(ng)+2,k)
              ad_A(-2,Mm(ng)+2,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+2,k)=tl_A(Lm(ng)-1, 2,k)
!>
              ad_A(Lm(ng)-1, 2,k)=ad_A(Lm(ng)-1, 2,k)+                  &
     &                            ad_A(-1,Mm(ng)+2,k)
              ad_A(-1,Mm(ng)+2,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+2,k)=tl_A(Lm(ng)  , 2,k)
!>
              ad_A(Lm(ng)  , 2,k)=ad_A(Lm(ng)  , 2,k)+                  &
     &                            ad_A( 0,Mm(ng)+2,k)
              ad_A( 0,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(-2,Mm(ng)+3,k)=tl_A(Lm(ng)-2, 3,k)
!>
                ad_A(Lm(ng)-2, 3,k)=ad_A(Lm(ng)-2, 3,k)+                &
     &                              ad_A(-2,Mm(ng)+3,k)
                ad_A(-2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(-1,Mm(ng)+3,k)=tl_A(Lm(ng)-1, 3,k)
!>
                ad_A(Lm(ng)-1, 3,k)=ad_A(Lm(ng)-1, 3,k)+                &
     &                              ad_A(-1,Mm(ng)+3,k)
                ad_A(-1,Mm(ng)+3,k)=0.0_r8
!>              tl_A( 0,Mm(ng)+3,k)=tl_A(Lm(ng)  , 3,k)
!>
                ad_A(Lm(ng)  , 3,k)=ad_A(Lm(ng)  , 3,k)+                &
     &                              ad_A( 0,Mm(ng)+3,k)
                ad_A( 0,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,Mm(ng)+1,k)=tl_A( 1, 1,k)
!>
              ad_A( 1, 1,k)=ad_A( 1, 1,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+1,k)
              ad_A(Lm(ng)+1,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+1,Mm(ng)+2,k)=tl_A( 1, 2,k)
!>
              ad_A( 1, 2,k)=ad_A( 1, 2,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+2,k)
              ad_A(Lm(ng)+1,Mm(ng)+2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+1,k)=tl_A( 2, 1,k)
!>
              ad_A( 2, 1,k)=ad_A( 2, 1,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+1,k)
              ad_A(Lm(ng)+2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+2,k)=tl_A( 2, 2,k)
!>
              ad_A( 2, 2,k)=ad_A( 2, 2,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+2,k)
              ad_A(Lm(ng)+2,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+1,Mm(ng)+3,k)=tl_A( 1, 3,k)
!>
                ad_A( 1, 3,k)=ad_A( 1, 3,k)+                            &
     &                        ad_A(Lm(ng)+1,Mm(ng)+3,k)
                ad_A(Lm(ng)+1,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+2,Mm(ng)+3,k)=tl_A( 2, 3,k)
!>
                ad_A( 2, 3,k)=ad_A( 2, 3,k)+                            &
     &                        ad_A(Lm(ng)+2,Mm(ng)+3,k)
                ad_A(Lm(ng)+2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+1,k)=tl_A( 3, 1,k)
!>
                ad_A( 3, 1,k)=ad_A( 3, 1,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+1,k)
                ad_A(Lm(ng)+3,Mm(ng)+1,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+2,k)=tl_A( 3, 2,k)
!>
                ad_A( 3, 2,k)=ad_A( 3, 2,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+2,k)
                ad_A(Lm(ng)+3,Mm(ng)+2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+3,k)=tl_A( 3, 3,k)
!>
                ad_A( 3, 3,k)=ad_A( 3, 3,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+3,k)
                ad_A(Lm(ng)+3,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
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
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,-2,k)=tl_A(i,Mm(ng)-2,k)
!>
                ad_A(i,Mm(ng)-2,k)=ad_A(i,Mm(ng)-2,k)+                  &
     &                             ad_A(i,-2,k)
                ad_A(i,-2,k)=0.0_r8
!>              tl_A(i,-1,k)=tl_A(i,Mm(ng)-1,k)
!>
                ad_A(i,Mm(ng)-1,k)=ad_A(i,Mm(ng)-1,k)+                  &
     &                             ad_A(i,-1,k)
                ad_A(i,-1,k)=0.0_r8
!>              tl_A(i, 0,k)=tl_A(i,Mm(ng)  ,k)
!>
                ad_A(i,Mm(ng)  ,k)=ad_A(i,Mm(ng)  ,k)+                  &
     &                             ad_A(i, 0,k)
                ad_A(i, 0,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,Mm(ng)+1,k)=tl_A(i,1,k)
!>
                ad_A(i,1,k)=ad_A(i,1,k)+                                &
     &                      ad_A(i,Mm(ng)+1,k)
                ad_A(i,Mm(ng)+1,k)=0.0_r8
!>              tl_A(i,Mm(ng)+2,k)=tl_A(i, 2,k)
!>
                ad_A(i,2,k)=ad_A(i,2,k)+                                &
     &                      ad_A(i,Mm(ng)+2,k)
                ad_A(i,Mm(ng)+2,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
!>                tl_A(i,Mm(ng)+3,k)=tl_A(i,3,k)
!>
                  ad_A(i,3,k)=ad_A(i,3,k)+                              &
     &                        ad_A(i,Mm(ng)+3,k)
                  ad_A(i,Mm(ng)+3,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
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
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(-2,j,k)=tl_A(Lm(ng)-2,j,k)
!>
                ad_A(Lm(ng)-2,j,k)=ad_A(Lm(ng)-2,j,k)+                  &
     &                             ad_A(-2,j,k)
                ad_A(-2,j,k)=0.0_r8
!>              tl_A(-1,j,k)=tl_A(Lm(ng)-1,j,k)
!>
                ad_A(Lm(ng)-1,j,k)=ad_A(Lm(ng)-1,j,k)+                  &
     &                             ad_A(-1,j,k)
                ad_A(-1,j,k)=0.0_r8
!>              tl_A( 0,j,k)=tl_A(Lm(ng)  ,j,k)
!>
                ad_A(Lm(ng)  ,j,k)=ad_A(Lm(ng)  ,j,k)+                  &
     &                             ad_A( 0,j,k)
                ad_A( 0,j,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(Lm(ng)+1,j,k)=tl_A(1,j,k)
!>
                ad_A(1,j,k)=ad_A(1,j,k)+                                &
     &                      ad_A(Lm(ng)+1,j,k)
                ad_A(Lm(ng)+1,j,k)=0.0_r8
!>              tl_A(Lm(ng)+2,j,k)=tl_A( 2,j,k)
!>
                ad_A(2,j,k)=ad_A(2,j,k)+                                &
     &                      ad_A(Lm(ng)+2,j,k)
                ad_A(Lm(ng)+2,j,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
!>                tl_A(Lm(ng)+3,j,k)=tl_A(3,j,k)
!>
                  ad_A(3,j,k)=ad_A(3,j,k)+                              &
     &                        ad_A(Lm(ng)+3,j,k)
                  ad_A(Lm(ng)+3,j,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_exchange_v3d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_w3d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
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
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,-2,k)=tl_A(Lm(ng)-2,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-2,k)=ad_A(Lm(ng)-2,Mm(ng)-2,k)+      &
     &                                  ad_A(-2,-2,k)
              ad_A(-2,-2,k)=0.0_r8
!>            tl_A(-2,-1,k)=tl_A(Lm(ng)-2,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)-1,k)=ad_A(Lm(ng)-2,Mm(ng)-1,k)+      &
     &                                  ad_A(-2,-1,k)
              ad_A(-2,-1,k)=0.0_r8
!>            tl_A(-2, 0,k)=tl_A(Lm(ng)-2,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-2,Mm(ng)  ,k)=ad_A(Lm(ng)-2,Mm(ng)  ,k)+      &
     &                                  ad_A(-2, 0,k)
              ad_A(-2, 0,k)=0.0_r8
!>            tl_A(-1,-2,k)=tl_A(Lm(ng)-1,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-2,k)=ad_A(Lm(ng)-1,Mm(ng)-2,k)+      &
     &                                  ad_A(-1,-2,k)
              ad_A(-1,-2,k)=0.0_r8
!>            tl_A(-1,-1,k)=tl_A(Lm(ng)-1,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)-1,k)=ad_A(Lm(ng)-1,Mm(ng)-1,k)+      &
     &                                  ad_A(-1,-1,k)
              ad_A(-1,-1,k)=0.0_r8
!>            tl_A(-1, 0,k)=tl_A(Lm(ng)-1,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)-1,Mm(ng)  ,k)=ad_A(Lm(ng)-1,Mm(ng)  ,k)+      &
     &                                  ad_A(-1, 0,k)
              ad_A(-1, 0,k)=0.0_r8
!>            tl_A( 0,-2,k)=tl_A(Lm(ng)  ,Mm(ng)-2,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-2,k)=ad_A(Lm(ng)  ,Mm(ng)-2,k)+      &
     &                                  ad_A( 0,-2,k)
              ad_A( 0,-2,k)=0.0_r8
!>            tl_A( 0,-1,k)=tl_A(Lm(ng)  ,Mm(ng)-1,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)-1,k)=ad_A(Lm(ng)  ,Mm(ng)-1,k)+      &
     &                                  ad_A( 0,-1,k)
              ad_A( 0,-1,k)=0.0_r8
!>            tl_A( 0, 0,k)=tl_A(Lm(ng)  ,Mm(ng)  ,k)
!>
              ad_A(Lm(ng)  ,Mm(ng)  ,k)=ad_A(Lm(ng)  ,Mm(ng)  ,k)+      &
     &                                  ad_A( 0, 0,k)
              ad_A( 0, 0,k)=0.0_r8
            END DO
          END IF

          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,-2,k)=tl_A( 1,Mm(ng)-2,k)
!>
              ad_A( 1,Mm(ng)-2,k)=ad_A( 1,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+1,-2,k)
              ad_A(Lm(ng)+1,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+1,-1,k)=tl_A( 1,Mm(ng)-1,k)
!>
              ad_A( 1,Mm(ng)-1,k)=ad_A( 1,Mm(ng)-1,k)+                  &
     &                          ad_A(Lm(ng)+1,-1,k)
              ad_A(Lm(ng)+1,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+1, 0,k)=tl_A( 1,Mm(ng)  ,k)
!>
              ad_A( 1,Mm(ng)  ,k)=ad_A( 1,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+1, 0,k)
              ad_A(Lm(ng)+1, 0,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-2,k)=tl_A( 2,Mm(ng)-2,k)
!>
              ad_A( 2,Mm(ng)-2,k)=ad_A( 2,Mm(ng)-2,k)+                  &
     &                            ad_A(Lm(ng)+2,-2,k)
              ad_A(Lm(ng)+2,-2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,-1,k)=tl_A( 2,Mm(ng)-1,k)
!>
              ad_A( 2,Mm(ng)-1,k)=ad_A( 2,Mm(ng)-1,k)+                  &
     &                            ad_A(Lm(ng)+2,-1,k)
              ad_A(Lm(ng)+2,-1,k)=0.0_r8
!>            tl_A(Lm(ng)+2, 0,k)=tl_A( 2,Mm(ng)  ,k)
!>
              ad_A( 2,Mm(ng)  ,k)=ad_A( 2,Mm(ng)  ,k)+                  &
     &                            ad_A(Lm(ng)+2, 0,k)
              ad_A(Lm(ng)+2, 0,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+3,-2,k)=tl_A(3 ,Mm(ng)-2,k)
!>
                ad_A(3 ,Mm(ng)-2,k)=ad_A(3 ,Mm(ng)-2,k)+                &
     &                              ad_A(Lm(ng)+3,-2,k)
                ad_A(Lm(ng)+3,-2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,-1,k)=tl_A(3 ,Mm(ng)-1,k)
!>
                ad_A(3 ,Mm(ng)-1,k)=ad_A(3 ,Mm(ng)-1,k)+                &
     &                              ad_A(Lm(ng)+3,-1,k)
                ad_A(Lm(ng)+3,-1,k)=0.0_r8
!>              tl_A(Lm(ng)+3, 0,k)=tl_A(3 ,Mm(ng)  ,k)
!>
                ad_A(3 ,Mm(ng)  ,k)=ad_A(3 ,Mm(ng)  ,k)+                &
     &                              ad_A(Lm(ng)+3, 0,k)
                ad_A(Lm(ng)+3, 0,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(-2,Mm(ng)+1,k)=tl_A(Lm(ng)-2, 1,k)
!>
              ad_A(Lm(ng)-2, 1,k)=ad_A(Lm(ng)-2, 1,k)+                  &
     &                            ad_A(-2,Mm(ng)+1,k)
              ad_A(-2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+1,k)=tl_A(Lm(ng)-1, 1,k)
!>
              ad_A(Lm(ng)-1, 1,k)=ad_A(Lm(ng)-1, 1,k)+                  &
     &                            ad_A(-1,Mm(ng)+1,k)
              ad_A(-1,Mm(ng)+1,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+1,k)=tl_A(Lm(ng)  , 1,k)
!>
              ad_A(Lm(ng)  , 1,k)=ad_A(Lm(ng)  , 1,k)+                  &
     &                            ad_A( 0,Mm(ng)+1,k)
              ad_A( 0,Mm(ng)+1,k)=0.0_r8
!>            tl_A(-2,Mm(ng)+2,k)=tl_A(Lm(ng)-2, 2,k)
!>
              ad_A(Lm(ng)-2, 2,k)=ad_A(Lm(ng)-2, 2,k)+                  &
     &                            ad_A(-2,Mm(ng)+2,k)
              ad_A(-2,Mm(ng)+2,k)=0.0_r8
!>            tl_A(-1,Mm(ng)+2,k)=tl_A(Lm(ng)-1, 2,k)
!>
              ad_A(Lm(ng)-1, 2,k)=ad_A(Lm(ng)-1, 2,k)+                  &
     &                            ad_A(-1,Mm(ng)+2,k)
              ad_A(-1,Mm(ng)+2,k)=0.0_r8
!>            tl_A( 0,Mm(ng)+2,k)=tl_A(Lm(ng)  , 2,k)
!>
              ad_A(Lm(ng)  , 2,k)=ad_A(Lm(ng)  , 2,k)+                  &
     &                            ad_A( 0,Mm(ng)+2,k)
              ad_A( 0,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(-2,Mm(ng)+3,k)=tl_A(Lm(ng)-2, 3,k)
!>
                ad_A(Lm(ng)-2, 3,k)=ad_A(Lm(ng)-2, 3,k)+                &
     &                              ad_A(-2,Mm(ng)+3,k)
                ad_A(-2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(-1,Mm(ng)+3,k)=tl_A(Lm(ng)-1, 3,k)
!>
                ad_A(Lm(ng)-1, 3,k)=ad_A(Lm(ng)-1, 3,k)+                &
     &                              ad_A(-1,Mm(ng)+3,k)
                ad_A(-1,Mm(ng)+3,k)=0.0_r8
!>              tl_A( 0,Mm(ng)+3,k)=tl_A(Lm(ng)  , 3,k)
!>
                ad_A(Lm(ng)  , 3,k)=ad_A(Lm(ng)  , 3,k)+                &
     &                              ad_A( 0,Mm(ng)+3,k)
                ad_A( 0,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
          END IF

          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
!>            tl_A(Lm(ng)+1,Mm(ng)+1,k)=tl_A( 1, 1,k)
!>
              ad_A( 1, 1,k)=ad_A( 1, 1,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+1,k)
              ad_A(Lm(ng)+1,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+1,Mm(ng)+2,k)=tl_A( 1, 2,k)
!>
              ad_A( 1, 2,k)=ad_A( 1, 2,k)+                              &
     &                      ad_A(Lm(ng)+1,Mm(ng)+2,k)
              ad_A(Lm(ng)+1,Mm(ng)+2,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+1,k)=tl_A( 2, 1,k)
!>
              ad_A( 2, 1,k)=ad_A( 2, 1,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+1,k)
              ad_A(Lm(ng)+2,Mm(ng)+1,k)=0.0_r8
!>            tl_A(Lm(ng)+2,Mm(ng)+2,k)=tl_A( 2, 2,k)
!>
              ad_A( 2, 2,k)=ad_A( 2, 2,k)+                              &
     &                      ad_A(Lm(ng)+2,Mm(ng)+2,k)
              ad_A(Lm(ng)+2,Mm(ng)+2,k)=0.0_r8
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
!>              tl_A(Lm(ng)+1,Mm(ng)+3,k)=tl_A( 1, 3,k)
!>
                ad_A( 1, 3,k)=ad_A( 1, 3,k)+                            &
     &                        ad_A(Lm(ng)+1,Mm(ng)+3,k)
                ad_A(Lm(ng)+1,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+2,Mm(ng)+3,k)=tl_A( 2, 3,k)
!>
                ad_A( 2, 3,k)=ad_A( 2, 3,k)+                            &
     &                        ad_A(Lm(ng)+2,Mm(ng)+3,k)
                ad_A(Lm(ng)+2,Mm(ng)+3,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+1,k)=tl_A( 3, 1,k)
!>
                ad_A( 3, 1,k)=ad_A( 3, 1,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+1,k)
                ad_A(Lm(ng)+3,Mm(ng)+1,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+2,k)=tl_A( 3, 2,k)
!>
                ad_A( 3, 2,k)=ad_A( 3, 2,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+2,k)
                ad_A(Lm(ng)+3,Mm(ng)+2,k)=0.0_r8
!>              tl_A(Lm(ng)+3,Mm(ng)+3,k)=tl_A( 3, 3,k)
!>
                ad_A( 3, 3,k)=ad_A( 3, 3,k)+                            &
     &                        ad_A(Lm(ng)+3,Mm(ng)+3,k)
                ad_A(Lm(ng)+3,Mm(ng)+3,k)=0.0_r8
              END DO
            END IF
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
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,-2,k)=tl_A(i,Mm(ng)-2,k)
!>
                ad_A(i,Mm(ng)-2,k)=ad_A(i,Mm(ng)-2,k)+                  &
     &                             ad_A(i,-2,k)
                ad_A(i,-2,k)=0.0_r8
!>              tl_A(i,-1,k)=tl_A(i,Mm(ng)-1,k)
!>
                ad_A(i,Mm(ng)-1,k)=ad_A(i,Mm(ng)-1,k)+                  &
     &                             ad_A(i,-1,k)
                ad_A(i,-1,k)=0.0_r8
!>              tl_A(i, 0,k)=tl_A(i,Mm(ng)  ,k)
!>
                ad_A(i,Mm(ng)  ,k)=ad_A(i,Mm(ng)  ,k)+                  &
     &                             ad_A(i, 0,k)
                ad_A(i, 0,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
!>              tl_A(i,Mm(ng)+1,k)=tl_A(i,1,k)
!>
                ad_A(i,1,k)=ad_A(i,1 ,k)+                               &
     &                      ad_A(i,Mm(ng)+1,k)
                ad_A(i,Mm(ng)+1,k)=0.0_r8
!>              tl_A(i,Mm(ng)+2,k)=tl_A(i,2,k)
!>
                ad_A(i,2,k)=ad_A(i,2,k)+                                &
     &                      ad_A(i,Mm(ng)+2,k)
                ad_A(i,Mm(ng)+2,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
!>                tl_A(i,Mm(ng)+3,k)=tl_A(i,3,k)
!>
                  ad_A(i,3,k)=ad_A(i,3,k)+                              &
     &                        ad_A(i,Mm(ng)+3,k)
                  ad_A(i,Mm(ng)+3,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
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
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(-2,j,k)=tl_A(Lm(ng)-2,j,k)
!>
                ad_A(Lm(ng)-2,j,k)=ad_A(Lm(ng)-2,j,k)+                  &
     &                             ad_A(-2,j,k)
                ad_A(-2,j,k)=0.0_r8
!>              tl_A(-1,j,k)=tl_A(Lm(ng)-1,j,k)
!>
                ad_A(Lm(ng)-1,j,k)=ad_A(Lm(ng)-1,j,k)+                  &
     &                             ad_A(-1,j,k)
                ad_A(-1,j,k)=0.0_r8
!>              tl_A( 0,j,k)=tl_A(Lm(ng)  ,j,k)
!>
                ad_A(Lm(ng)  ,j,k)=ad_A(Lm(ng)  ,j,k)+                  &
     &                             ad_A( 0,j,k)
                ad_A( 0,j,k)=0.0_r8
              END DO
            END DO
          END IF

          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
!>              tl_A(Lm(ng)+1,j,k)=tl_A(1,j,k)
!>
                ad_A(1,j,k)=ad_A(1,j,k)+                                &
     &                      ad_A(Lm(ng)+1,j,k)
                ad_A(Lm(ng)+1,j,k)=0.0_r8
!>              tl_A(Lm(ng)+2,j,k)=tl_A( 2,j,k)
!>
                ad_A(2,j,k)=ad_A(2,j,k)+                                &
     &                      ad_A(Lm(ng)+2,j,k)
                ad_A(Lm(ng)+2,j,k)=0.0_r8
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
!>                tl_A(Lm(ng)+3,j,k)=tl_A(3,j,k)
!>
                  ad_A(3,j,k)=ad_A(3,j,k)+                              &
     &                        ad_A(Lm(ng)+3,j,k)
                  ad_A(Lm(ng)+3,j,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_exchange_w3d_tile
#endif
      END MODULE ad_exchange_3d_mod
