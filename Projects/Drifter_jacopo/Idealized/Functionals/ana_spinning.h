      SUBROUTINE ana_spinning (ng, tile, model)
!
!! svn $Id: ana_spinning.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This subroutine sets time-variable rotation force as the sum of     !
!  Coriolis and Centripetal accelerations.  This is used in polar      !
!  coordinate applications (annulus grid).                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_spinning_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
#ifdef SPHERICAL
     &                        GRID(ng) % lonr,                          &
     &                        GRID(ng) % latr,                          &
#else
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr,                            &
#endif
     &                        GRID(ng) % f,                             &
     &                        GRID(ng) % omn,                           &
     &                        GRID(ng) % fomn)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(26),'(a,a)') TRIM(Adir), '/ana_spinning.h'
      END IF

      RETURN
      END SUBROUTINE ana_spinning
!
!***********************************************************************
      SUBROUTINE ana_spinning_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
#ifdef SPHERICAL
     &                              lonr, latr                          &
#else
     &                              xr, yr,                             &
#endif
     &                              f, omn, fomn)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: fomn(LBi:,LBj:)
#else
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: fomn(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time-varying rotation force: Coriolis plus Centripetal
!  accelerations.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          fomn(i,j)=???
        END DO
      END DO
#else
      ana_spinningr: No values provided for Coriolis + Centripetal
                     accelerations.
#endif

      RETURN
      END SUBROUTINE ana_spinning_tile
