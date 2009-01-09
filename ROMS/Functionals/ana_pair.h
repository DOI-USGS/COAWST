      SUBROUTINE ana_pair (ng, tile, model)
!
!! svn $Id: ana_pair.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
     &                    FORCES(ng) % Pair)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(17)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_pair
!
!***********************************************************************
      SUBROUTINE ana_pair_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Pair)
!***********************************************************************
!
      USE mod_param
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Pair(LBi:,LBj:)
#else
      real(r8), intent(out) :: Pair(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air pressure (mb).
!  (1 mb = 100 Pa = 1 hPa,  1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).
!-----------------------------------------------------------------------
!
#if defined BENCHMARK
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Pair(i,j)=1025.0_r8
        END DO
      END DO
#elif defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Pair(i,j)=1013.48_r8
        END DO
      END DO
#else
      ana_pair.h: no values provided for Pair.
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pair)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Pair)
#endif
      RETURN
      END SUBROUTINE ana_pair_tile
