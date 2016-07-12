#include "cppdefs.h"
      MODULE ad_bc_2d_mod
#ifdef ADJOINT
!
!svn $Id: ad_bc_2d.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 2D adjoint fields.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     ad_A              2D adjoint field.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A              Processed 2D adjoint field.                    !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     ad_bc_r2d_tile    Boundary conditions for field at RHO-points    !
!     ad_bc_u2d_tile    Boundary conditions for field at U-points      !
!     ad_bc_v2d_tile    Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bc_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+            &
!>   &                                  tl_A(Istr  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jend+1)
            ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
            ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
            ad_A(Istr-1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+            &
!>   &                                  tl_A(Iend+1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
            ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
            ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
            ad_A(Iend+1,Jstr-1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          tl_A(Istr-1,Jstr-1)=0.5_r8*(tl_A(Istr  ,Jstr-1)+            &
!>   &                                  tl_A(Istr-1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jstr-1)
            ad_A(Istr  ,Jstr-1)=ad_A(Istr  ,Jstr-1)+adfac
            ad_A(Istr-1,Jstr  )=ad_A(Istr-1,Jstr  )+adfac
            ad_A(Istr-1,Jstr-1)=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
              ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
              ad_A(i,Jstr-1)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            tl_A(i,Jend+1)=tl_A(i,Jend)
!>
              ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
              ad_A(i,Jend+1)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            tl_A(Istr-1,j)=tl_A(Istr,j)
!>
              ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
              ad_A(Istr-1,j)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            tl_A(Iend+1,j)=tl_A(Iend,j)
!>
              ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
              ad_A(Iend+1,j)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_u2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr  ,Jend+1)=0.5_r8*(tl_A(Istr  ,Jend  )+            &
!>   &                                  tl_A(Istr+1,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr,Jend+1)
            ad_A(Istr  ,Jend  )=ad_A(Istr  ,Jend  )+adfac
            ad_A(Istr+1,Jend+1)=ad_A(Istr+1,Jend+1)+adfac
            ad_A(Istr  ,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+            &
!>   &                                  tl_A(Iend+1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
            ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
            ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
            ad_A(Iend+1,Jstr-1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          tl_A(Istr  ,Jstr-1)=0.5_r8*(tl_A(Istr+1,Jstr-1)+            &
!>   &                                  tl_A(Istr  ,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Istr,Jstr-1)
            ad_A(Istr+1,Jstr-1)=ad_A(Istr+1,Jstr-1)+adfac
            ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
            ad_A(Istr  ,Jstr-1)=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed (free-slip/no-slip)
!  or gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (ad_LBC(isouth,isBu2d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO i=Imin,Imax
              IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>              tl_A(i,Jstr-1)=tl_A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
!>
                ad_A(i,Jstr-1)=ad_A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
# endif
!>              tl_A(i,Jstr-1)=gamma2(ng)*tl_A(i,Jstr)
!>
                ad_A(i,Jstr)=ad_A(i,Jstr)+gamma2(ng)*ad_A(i,Jstr-1)
                ad_A(i,Jstr-1)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
                ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
                ad_A(i,Jstr-1)=0.0_r8
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (ad_LBC(inorth,isBu2d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO i=Imin,Imax
              IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>              tl_A(i,Jend+1)=tl_A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
!>
                ad_A(i,Jend+1)=ad_A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
# endif
!>              tl_A(i,Jend+1)=gamma2(ng)*tl_A(i,Jend)
!>
                ad_A(i,Jend)=ad_A(i,Jend)+gamma2(ng)*ad_A(i,Jend+1)
                ad_A(i,Jend+1)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              tl_A(i,Jend+1)=tl_A(i,Jend)
!>
                ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
                ad_A(i,Jend+1)=0.0_r8
              END IF
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (ad_LBC(iwest,isBu2d,ng)%closed) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              tl_A(Istr,j)=0.0_r8
!>
                ad_A(Istr,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              tl_A(Istr,j)=tl_A(Istr+1,j)
!>
                ad_A(Istr+1,j)=ad_A(Istr+1,j)+ad_A(Istr,j)
                ad_A(Istr  ,j)=0.0_r8
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (ad_LBC(ieast,isBu2d,ng)%closed) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              A(Iend+1,j)=0.0_r8
!>
                ad_A(Iend+1,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              tl_A(Iend+1,j)=tl_A(Iend,j)
!>
                ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
                ad_A(Iend+1,j)=0.0_r8
              END IF
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_v2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+            &
!>   &                                  tl_A(Istr  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jend+1)
            ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
            ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
            ad_A(Istr-1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
!>          tl_A(Iend+1,Jstr  )=0.5_r8*(tl_A(Iend  ,Jstr  )+            &
!>   &                                  tl_A(Iend+1,Jstr+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr  )
            ad_A(Iend  ,Jstr  )=ad_A(Iend  ,Jstr  )+adfac
            ad_A(Iend+1,Jstr+1)=ad_A(Iend+1,Jstr+1)+adfac
            ad_A(Iend+1,Jstr  )=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
!>          tl_A(Istr-1,Jstr  )=0.5_r8*(tl_A(Istr  ,Jstr  )+            &
!>   &                                  tl_A(Istr-1,Jstr+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jstr  )
            ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
            ad_A(Istr-1,Jstr+1)=ad_A(Istr-1,Jstr+1)+adfac
            ad_A(Istr-1,Jstr  )=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (ad_LBC(isouth,isBv2d,ng)%closed) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              tl_A(i,Jstr)=0.0_r8
!>
                ad_A(i,Jstr)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              tl_A(i,Jstr)=tl_A(i,Jstr+1)
!>
                ad_A(i,Jstr+1)=ad_A(i,Jstr+1)+ad_A(i,Jstr)
                ad_A(i,Jstr)=0.0_r8
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (ad_LBC(inorth,isBv2d,ng)%closed) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              tl_A(i,Jend+1)=0.0_r8
!>
                ad_A(i,Jend+1)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              tl_A(i,Jend+1)=tl_A(i,Jend)
!>
                ad_A(i,Jend)=ad_A(i,Jend)+ad_A(i,Jend+1)
                ad_A(i,Jend+1)=0.0_r8
              END IF
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (ad_LBC(iwest,isBv2d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>              tl_A(Istr-1,j)=tl_A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
!>
                ad_A(Istr-1,j)=ad_A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
# endif
!>              tl_A(Istr-1,j)=gamma2(ng)*tl_A(Istr,j)
!>
                ad_A(Istr  ,j)=ad_A(Istr,j)+gamma2(ng)*ad_A(Istr-1,j)
                ad_A(Istr-1,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              tl_A(Istr-1,j)=tl_A(Istr,j)
!>
                ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
                ad_A(Istr-1,j)=0.0_r8
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (ad_LBC(ieast,isBv2d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>              tl_A(Iend+1,j)=tl_A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
!>
                ad_A(Iend+1,j)=ad_A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
# endif
!>              tl_A(Iend+1,j)=gamma2(ng)*tl_A(Iend,j)
!>
                ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
                ad_A(Iend+1,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              tl_A(Iend+1,j)=tl_A(Iend,j)
!>
                ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
                ad_A(Iend+1,j)=0.0_r8
              END IF
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_v2d_tile

!
!***********************************************************************
      SUBROUTINE ad_dabc_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+            &
!>   &                                  tl_A(Istr  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jend+1)
            ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
            ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
            ad_A(Istr-1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+            &
!>   &                                  tl_A(Iend+1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
            ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
            ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
            ad_A(Iend+1,Jstr-1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          tl_A(Istr-1,Jstr-1)=0.5_r8*(tl_A(Istr  ,Jstr-1)+            &
!>   &                                  tl_A(Istr-1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jstr-1)
            ad_A(Istr  ,Jstr-1)=ad_A(Istr  ,Jstr-1)+adfac
            ad_A(Istr-1,Jstr  )=ad_A(Istr-1,Jstr  )+adfac
            ad_A(Istr-1,Jstr-1)=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
              ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
              ad_A(i,Jstr-1)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            tl_A(i,Jend+1)=tl_A(i,Jend)
!>
              ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
              ad_A(i,Jend+1)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            tl_A(Istr-1,j)=tl_A(Istr,j)
!>
              ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
              ad_A(Istr-1,j)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            tl_A(Iend+1,j)=tl_A(Iend,j)
!>
              ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
              ad_A(Iend+1,j)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_dabc_r2d_tile

!
!***********************************************************************
      SUBROUTINE ad_dabc_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_u2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr  ,Jend+1)=0.5_r8*(tl_A(Istr  ,Jend  )+            &
!>   &                                  tl_A(Istr+1,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr,Jend+1)
            ad_A(Istr  ,Jend  )=ad_A(Istr  ,Jend  )+adfac
            ad_A(Istr+1,Jend+1)=ad_A(Istr+1,Jend+1)+adfac
            ad_A(Istr  ,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+            &
!>   &                                  tl_A(Iend+1,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
            ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
            ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
            ad_A(Iend+1,Jstr-1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          tl_A(Istr  ,Jstr-1)=0.5_r8*(tl_A(Istr+1,Jstr-1)+            &
!>   &                                  tl_A(Istr  ,Jstr  ))
!>
            adfac=0.5_r8*ad_A(Istr,Jstr-1)
            ad_A(Istr+1,Jstr-1)=ad_A(Istr+1,Jstr-1)+adfac
            ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
            ad_A(Istr  ,Jstr-1)=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
              ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
              ad_A(i,Jstr-1)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            tl_A(i,Jend+1)=tl_A(i,Jend)
!>
              ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
              ad_A(i,Jend+1)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            tl_A(Istr,j)=tl_A(Istr+1,j)
!>
              ad_A(Istr+1,j)=ad_A(Istr+1,j)+ad_A(Istr,j)
              ad_A(Istr  ,j)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            tl_A(Iend+1,j)=tl_A(Iend,j)
!>
              ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
              ad_A(Iend+1,j)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_dabc_u2d_tile

!
!***********************************************************************
      SUBROUTINE ad_dabc_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_v2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+            &
!>   &                                  tl_A(Iend  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jend+1)
            ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
            ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
            ad_A(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+            &
!>   &                                  tl_A(Istr  ,Jend+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jend+1)
            ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
            ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
            ad_A(Istr-1,Jend+1)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
!>          tl_A(Iend+1,Jstr  )=0.5_r8*(tl_A(Iend  ,Jstr  )+            &
!>   &                                  tl_A(Iend+1,Jstr+1))
!>
            adfac=0.5_r8*ad_A(Iend+1,Jstr  )
            ad_A(Iend  ,Jstr  )=ad_A(Iend  ,Jstr  )+adfac
            ad_A(Iend+1,Jstr+1)=ad_A(Iend+1,Jstr+1)+adfac
            ad_A(Iend+1,Jstr  )=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
!>          tl_A(Istr-1,Jstr  )=0.5_r8*(tl_A(Istr  ,Jstr  )+            &
!>   &                                  tl_A(Istr-1,Jstr+1))
!>
            adfac=0.5_r8*ad_A(Istr-1,Jstr  )
            ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
            ad_A(Istr-1,Jstr+1)=ad_A(Istr-1,Jstr+1)+adfac
            ad_A(Istr-1,Jstr  )=0.0_r8
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            tl_A(i,Jstr)=tl_A(i,Jstr+1)
!>
              ad_A(i,Jstr+1)=ad_A(i,Jstr+1)+ad_A(i,Jstr)
              ad_A(i,Jstr)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            tl_A(i,Jend+1)=tl_A(i,Jend)
!>
              ad_A(i,Jend)=ad_A(i,Jend)+ad_A(i,Jend+1)
              ad_A(i,Jend+1)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            tl_A(Istr-1,j)=tl_A(Istr,j)
!>
              ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
              ad_A(Istr-1,j)=0.0_r8
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            tl_A(Iend+1,j)=tl_A(Iend,j)
!>
              ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
              ad_A(Iend+1,j)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_dabc_v2d_tile
#endif
      END MODULE ad_bc_2d_mod
