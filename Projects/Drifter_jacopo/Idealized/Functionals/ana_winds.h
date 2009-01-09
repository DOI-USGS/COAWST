      SUBROUTINE ana_winds (ng, tile, model)
!
!! svn $Id: ana_winds.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
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
      CALL ana_winds_tile (ng, model, Istr, Iend, Jstr, Jend,           &
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
        WRITE (ANANAME(36),'(a,a)')  TRIM(Adir), 'ana_winds.h'
      END IF

      RETURN
      END SUBROUTINE ana_winds
!
!***********************************************************************
      SUBROUTINE ana_winds_tile (ng, model, Istr, Iend, Jstr, Jend,     &
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
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set surface wind components (m/s) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Uwind(i,j)=???
          Vwind(i,j)=???
        END DO
      END DO
#else
      ana_winds.h: no values provided for Uwind and Vwind.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Uwind)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Vwind)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Uwind, Vwind)
#endif

      RETURN
      END SUBROUTINE ana_winds_tile
