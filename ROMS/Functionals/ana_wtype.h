      SUBROUTINE ana_wtype (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets spatially varying Jerlov water type index.     !
!  It is used in 'lmd_swfrac' to compute the fraction of shortwave     !
!  flux penetrating the water column (light absorption), modeled       !
!  as a double exponential decay function in Jerlov water type.        !
!                                                                      !
!  Currently, the following Jerlov water types are supported:          !
!                                                                      !
!  Array     Jerlov                                                    !
!  Index   Water Type   Examples                                       !
!  -----   ----------   --------                                       !
!                                                                      !
!    1         I        Open Pacific                                   !
!    2         IA       Eastern Mediterranean, Indian Ocean            !
!    3         IB       Western Mediterranean, Open Atlantic           !
!    4         II       Coastal waters, Azores                         !
!    5         III      Coastal waters, North Sea                      !
!    6         1        Skagerrak Strait                               !
!    7         3        Baltic                                         !
!    8         5        Black Sea                                      !
!    9         7        Coastal waters, dark                           !
!                                                                      !
!  The range of indices 1:9 are ordered by increasing absorption:      !
!  from clear water (type I) to dark turbidity water (type 7).         !
!  The indices correspond to the paramenters used to model the         !
!  the light absorption into the water column using a double           !
!  exponential fitting function of (Paulson and Simpson, 1997).        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_wtype_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     GRID(ng) % h,                                &
     &                     MIXING(ng) % Jwtype)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(39)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_wtype
!
!***********************************************************************
      SUBROUTINE ana_wtype_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           h, Jwtype)
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
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(out) :: Jwtype(LBi:,LBj:)
#else
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Jwtype(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: fac

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Jerlov water type array indices (1 to 9, currently) for light
!  absorption.
!-----------------------------------------------------------------------
!
#ifdef MY_APPLICATION
      fac=1.0/1000.0_r8                   ! Inverse bathymetry threshold
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Jwtype(i,j)=ANINT(5.0_r8-4.5_r8*(TANH(h(i,j)*fac)))      ! 1:5
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Jwtype(i,j)=???
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Jwtype)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Jwtype)
#endif

      RETURN
      END SUBROUTINE ana_wtype_tile
