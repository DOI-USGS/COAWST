#include "cppdefs.h"
      MODULE wec_dissip_mod
#if defined SOLVE3D && (defined WDISS_THORGUZA || \
                        defined WDISS_CHURTHOR)
!svn $Id: wec_dissip.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the terms corresponding to vortex forces in   !
!  momentum equations.                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Thornton, E. B., and R. T. Guza, Surf zone longshore currents and   !
!  random waves: field data and models, J. Phys. Oceanogr.,            !
!  16,1165–1178, 1986.                                                 !
!                                                                      !
!  Church, J. C., and E. B. Thornton, Effects of breaking wave induced !
!  turbulence within a longshore current model, Coastal Eng., 20,      !
!  1–28, 1993.                                                         !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_dissip
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_dissip (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_param
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
!
      integer, intent(in) :: ng, tile
# include "tile.h"
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL wec_dissip_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),        &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          GRID(ng) % h,                           &
     &                          OCEAN(ng) % zeta,                       &
     &                          FORCES(ng) % Hwave,                     &
     &                          FORCES(ng) % Pwave_top,                 &
     &                          FORCES(ng) % Dissip_break,              &
     &                          FORCES(ng) % Dissip_wcap)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
      RETURN
      END SUBROUTINE wec_dissip
!
!***********************************************************************
      SUBROUTINE wec_dissip_tile (ng, tile,    LBi, UBi, LBj, UBj, UBk, &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            h, zeta,                              &
     &                            Hwave, Pwave_top,                     &
     &                            Dissip_break, Dissip_wcap)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcap(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcap(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff1, cff2
      real(r8) :: fac1, sigmat, Dstp
# if defined WDISS_CHURTHOR
      real(r8) :: RB1, RB2
# endif
      real(r8), parameter :: gammaw=0.31_r8
      real(r8), parameter :: eps = 1.0E-14_r8
# include "set_bounds.h"
!
      fac1=3.0_r8*g*sqrt(pi)/16.0_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Compute total depth
!
          Dstp=zeta(i,j,1)+h(i,j)
          cff1=0.707_r8*Hwave(i,j)
          sigmat=MIN(1.0_r8/(Pwave_top(i,j)+eps),1.0_r8)
!
# ifdef WDISS_THORGUZA
!
!  Calcualate wave dissipation using empirical parameters of 
!  Thornton and Guza, 1986
!
          cff2=1.0_r8/((gammaw**4.0_r8)*(Dstp**5.0_r8))
          Dissip_break(i,j)=0.2621_r8*fac1*sigmat*                      &
     &                     (cff1**7.0_r8)*cff2
          Dissip_wcap(i,j)=0.0_r8
!
# elif defined WDISS_CHURTHOR
!
!  Calculate wave dissipation using empirical parameters of
!  Church and Thornton, 1993.
! 
          cff2=1.0_r8/(gammaw*Dstp)
          RB1=1.0_r8+tanh(8.0_r8*((cff1*cff2)-1.0_r8))
          RB2=1.0_r8-(1.0_r8+(cff1*cff2)**2.0_r8)**(-2.5_r8) 
          Dissip_break(i,j)=(0.2621_r8/Dstp)*fac1*sigmat*               &
     &                      (cff1**3.0_r8)*RB1*RB2
          Dissip_wcap(i,j)=0.0_r8
# endif
        END DO
      END DO
# if defined WDISS_THORGUZA || defined WDISS_CHURTHOR
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_break)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcap)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break, Dissip_wcap)
# endif
      RETURN
      END SUBROUTINE wec_dissip_tile
#endif
      END MODULE wec_dissip_mod
