      SUBROUTINE ana_rain (ng, tile, model)
!
!! svn $Id: ana_rain.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets precipitation rate (kg/m2/s) using an             !
!  analytical expression.                                              !
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
      CALL ana_rain_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    FORCES(ng) % rain)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(21),'(a,a)') TRIM(Adir), '/ana_rain.h'
      END IF

      RETURN
      END SUBROUTINE ana_rain
!
!***********************************************************************
      SUBROUTINE ana_rain_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rain)
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
      real(r8), intent(out) :: rain(LBi:,LBj:)
#else
      real(r8), intent(out) :: rain(LBi:UBi,LBj:UBj)
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
!  Set analytical precipitation rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rain(i,j)=???
        END DO
      END DO
#else
      ana_rain.h: No values provided for rain.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rain)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rain)
#endif

      RETURN
      END SUBROUTINE ana_rain_tile
