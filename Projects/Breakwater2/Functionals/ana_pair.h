      SUBROUTINE ana_pair (ng, tile, model)
!
!! svn $Id: ana_pair.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
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
      CALL ana_pair_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    FORCES(ng) % Pair)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(17),'(a,a)') TRIM(Adir), '/ana_pair.h'
      END IF

      RETURN
      END SUBROUTINE ana_pair
!
!***********************************************************************
      SUBROUTINE ana_pair_tile (ng, model, Istr, Iend, Jstr, Jend,      &
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
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air pressure (mb).
!  (1 mb = 100 Pa = 1 hPa,  1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Pair(i,j)=???
        END DO
      END DO
#else
      ana_pair.h: No values provided for Pair.      
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pair)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Pair)
#endif

      RETURN
      END SUBROUTINE ana_pair_tile
