      SUBROUTINE ana_smflux (ng, tile, model)
!
!! svn $Id: ana_smflux.h 798 2007-05-21 19:44:03Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
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
      CALL ana_smflux_tile (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % angler,                          &
#ifdef SPHERICAL
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#else
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
#endif
#ifdef TL_IOMS
     &                      FORCES(ng) % tl_sustr,                      &
     &                      FORCES(ng) % tl_svstr,                      &
#endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(24)='ROMS/Functionals/ana_smflux.h'
      END IF

      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            angler,                               &
#ifdef SPHERICAL
     &                            lonr, latr,                           &
#else
     &                            xr, yr,                               &
#endif
#ifdef TL_IOMS
     &                            tl_sustr, tl_svstr,                   &
#endif
     &                            sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
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
      real(r8), intent(in) :: angler(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(out) :: tl_svstr(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_svstr(LBi:UBi,LBj:UBj)
# endif
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
      integer :: i, j, idum
      real(r8) :: Ewind, Nwind, cff, val1, val2, ran

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
#ifdef EKMAN_WALL
      DO j=JstrR,JendR 
        DO i=Istr,IendR 
          sustr(i,j)=0.0001_r8/rho0
        END DO 
      END DO 
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  ETA-direction (m2/s2) at horizontal V-points.
!-----------------------------------------------------------------------
!
#if defined EKMAN_WALL
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          svstr(i,j)=0.0_r8
        END DO
      END DO
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr)
# ifdef TL_IOMS
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_sustr)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_svstr)
# endif
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sustr, svstr)
#  ifdef TL_IOMS
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_sustr, tl_svstr)
#  endif
#endif
      RETURN
      END SUBROUTINE ana_smflux_tile
