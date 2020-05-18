      SUBROUTINE ana_pair (ng, tile, model)
!
!! svn $Id: ana_pair.h 995 2020-01-10 04:01:28Z arango $
!!======================================================================
!! Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets surface air pressure (mb) using an analytical     !
!  expression.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_pair_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    FORCES(ng) % Pair)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(17)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_pair
!
!***********************************************************************
      SUBROUTINE ana_pair_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Pair)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Pair(LBi:,LBj:)
#else
      real(r8), intent(out) :: Pair(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air pressure (mb).
!  (1 mb = 100 Pa = 1 hPa,  1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Pair(i,j)=???
        END DO
      END DO
#else
      ana_pair.h: no values provided for Pair.
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Pair)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Pair)
#endif

      RETURN
      END SUBROUTINE ana_pair_tile
