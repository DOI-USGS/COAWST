      SUBROUTINE ana_humid (ng, tile, model)
!
!! svn $Id: ana_humid.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets surface air humidity (moisture) using an          !
!  analytical expression.  There three types of humidity:              !
!                                                                      !
!     1) Absolute humidity: density of water vapor.                    !
!     2) Specific humidity: ratio of the mass of water vapor to        !
!        the mass of moist air cointaining the vapor (g/kg)            !
!     3) Relative humidity: ratio of the actual mixing ratio to        !
!        saturation mixing ratio of the air at given temperature       !
!        and pressure (percentage).                                    !
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
      CALL ana_humid_tile (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng) % Hair)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME( 9),'(a,a)') TRIM(Adir), '/ana_humid.h'
      END IF

      RETURN
      END SUBROUTINE ana_humid
!
!***********************************************************************
      SUBROUTINE ana_humid_tile (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Hair)
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
      real(r8), intent(out) :: Hair(LBi:,LBj:)
#else
      real(r8), intent(out) :: Hair(LBi:UBi,LBj:UBj)
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
!  Set analytical surface air humidity.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Hair(i,j)=???
        END DO
      END DO
#else
      ana_humidity.h: no values provided for Hair.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Hair)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Hair)
#endif

      RETURN
      END SUBROUTINE ana_humid_tile
