      SUBROUTINE ana_winds (ng, tile, model)
!
!! svn $Id: ana_winds.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets surface wind components using an analytical       !
!  expression.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_winds_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
#ifdef SPHERICAL
     &                     GRID(ng) % lonr,                             &
     &                     GRID(ng) % latr,                             &
#else
     &                     GRID(ng) % xr,                               &
     &                     GRID(ng) % yr,                               &
#endif
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(36)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_winds
!
!***********************************************************************
      SUBROUTINE ana_winds_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
#ifdef SPHERICAL
     &                           lonr, latr,                            &
#else
     &                           xr, yr,                                &
#endif
     &                           Uwind, Vwind)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
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
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: Uwind(LBi:,LBj:)
      real(r8), intent(out) :: Vwind(LBi:,LBj:)
#else
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Vwind(LBi:UBi,LBj:UBj)
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
      real(r8) :: Wdir, Wmag, cff, u_wind, v_wind

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set surface wind components (m/s) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined BENCHMARK
      Wmag=15.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=0.2_r8*(60.0_r8+latr(i,j))
          Uwind(i,j)=Wmag*EXP(-cff*cff)
          Vwind(i,j)=0.0_r8
        END DO
      END DO
#elif defined BL_TEST
      IF ((tdays(ng)-dstart).le.6.0_r8) THEN
        u_wind=0.0_r8
!!      v_wind=4.7936_r8
        v_wind=10.0_r8
      END IF
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Uwind(i,j)=u_wind
          Vwind(i,j)=v_wind
        END DO
      END DO
#else
      ana_winds.h: No values provided for Uwind and Vwind.
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Uwind)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Vwind)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Uwind, Vwind)
#endif
      RETURN
      END SUBROUTINE ana_winds_tile
