      SUBROUTINE ana_spinning (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
      CALL ana_spinning_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
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
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(26)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_spinning
!
!***********************************************************************
      SUBROUTINE ana_spinning_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
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
#ifdef LAB_CANYON
      real(r8), parameter :: Omega0 = 2.0_r8*pi/25.0_r8
      real(r8), parameter :: Width = 0.20_r8
      real(r8), parameter :: Ro = 0.10_r8
      real(r8), parameter :: Rs = 0.55_r8
      real(r8), parameter :: little_omega = 2.0_r8*pi/24.0_r8
      real(r8), parameter :: Bu = 10.0_r8
      real(r8), parameter :: hd = 0.125_r8

      real(r8) :: Omega1, Omega1_of_t, Ro_t
      real(r8) :: fcor, d_rho_dz, d_Omega1_dt, time_fac
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time-varying rotation force: Coriolis plus Centripetal
!  accelerations.
!-----------------------------------------------------------------------
!
#ifdef LAB_CANYON
      fcor=2.0_r8*Omega0
      Omega1=fcor*Width*Ro/Rs
      Ro_t=little_omega/fcor
      d_rho_dz=(1000.0_r8*Bu/g)*(fcor*Width/hd)**2
      time_fac=1.0_r8+(Omega1/Omega0)*SIN(little_omega*time(ng))
      Omega1_of_t=Omega1*SIN(little_omega*time(ng))
      d_Omega1_dt=Omega1*little_omega*COS(little_omega*time(ng))
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          fomn(i,j)=(f(i,j)*time_fac+                                   &
     &               SQRT(xr(i,j)*xr(i,j)+yr(i,j)*yr(i,j))*             &
     &               ((2.0_r8*Omega0+Omega1_of_t)*Omega1_of_t))*        &
     &              omn(i,j)
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE ana_spinning_tile
