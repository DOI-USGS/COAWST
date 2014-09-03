#include "cppdefs.h"
      MODULE wec_roller_mod
#if defined SOLVE3D && defined WEC_ROLLER
!
!svn $Id: wec_roller.F 1428 2008-03-12 13:07:21Z jcwarner $
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
!  Svendsen, I.A., 1984. Mass flux and undertow in a surf zone.        !
!  Coastal Engineering 8, pp. 347–365.                              !
!                                                                      !
!  Reniers, A.J.M.H., Roelvink, J.A., and Thornton, E.B., 2004.        !
!  Morphodynamic modeling of an embayed beach under wave group forcing.!
!   J. Geophys. Res., 109: C01030, doi:10.1029/2002JC001586.           !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_roller
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_roller (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
!
      integer, intent(in) :: ng, tile
# include "tile.h"
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL wec_roller_tile (ng, tile,    LBi, UBi, LBj, UBj, N(ng),     &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs(ng),                             &
     &                            GRID(ng) % angler,                    &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % Hz,                        &
# ifdef ROLLER_RENIERS
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % pm,                        &
     &                            GRID(ng) % pn,                        &
# endif
     &                            OCEAN(ng) % ubar,                     &
     &                            OCEAN(ng) % vbar,                     &
     &                            OCEAN(ng) % zeta,                     &
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
     &                            FORCES(ng) % Dissip_break,            &
# if defined ROLLER_SVENDSEN || ROLLER_MONO
     &                            FORCES(ng) % Wave_break,              &
# endif
     &                            FORCES(ng) % Dissip_roller,           &
     &                            FORCES(ng) % rollA)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
      RETURN
      END SUBROUTINE wec_roller
!
!***********************************************************************
      SUBROUTINE wec_roller_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,    &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  nrhs,                           &
     &                                  angler, h, Hz,                  &
# ifdef ROLLER_RENIERS
     &                                  on_u, om_v, pm, pn,             &
# endif
     &                                  ubar, vbar, zeta,               &
     &                                  Hwave, Dwave, Lwave,            &
     &                                  Dissip_break,                   &
# if defined ROLLER_SVENDSEN || ROLLER_MONO
     &                                  Wave_break,                     &
# endif                                   
     &                                  Dissip_roller,                  &
     &                                  rollA)
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
      integer, intent(in) :: nrhs
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#  ifdef ROLLER_RENIERS
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_break(LBi:,LBj:)
#  if defined ROLLER_SVENDSEN || ROLLER_MONO
      real(r8), intent(in) :: Wave_break(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Dissip_roller(LBi:,LBj:)
      real(r8), intent(inout) :: rollA(LBi:,LBj:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
#  ifdef ROLLER_RENIERS
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dissip_break(LBi:UBi,LBj:UBj)
#  if defined ROLLER_SVENDSEN || ROLLER_MONO
      real(r8), intent(in) :: Wave_break(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Dissip_roller(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rollA(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, numits, it
      real(r8) :: cff, cff1, cff2, cff3
      real(r8), parameter :: sinb=0.1_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
# include "set_bounds.h"
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=zeta(i,j,1)+h(i,j)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j)))
          osigma(i,j)=1.0_r8/sigma(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
          wavec(i,j)=SQRT(g*owaven(i,j)*TANH(kD(i,j)))          
        END DO
      END DO
!
# if defined ROLLER_SVENDSEN
!
!     Check if the calculation of Dissip Roller is Consistent
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff1=0.0424_r8*Hwave(i,j)*Wave_break(i,j)
          rollA(i,j)=cff1*wavec(i,j)*wavec(i,j)*osigma(i,j)
          Dissip_roller(i,j)=g*sinb*rollA(i,j)*sigma(i,j)/wavec(i,j)
        END DO
      END DO
# elif defined ROLLER_MONO
!
!     Here Wave_break is really Breaking Area.
!     Check if the calculation of Dissip Roller is Consistent
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff1=Wave_break(i,j)/MAX(Lwave(i,j),Lwave_min)
          rollA(i,j)=cff1*wavec(i,j)*wavec(i,j)*osigma(i,j)
          Dissip_roller(i,j)=g*sinb*rollA(i,j)*sigma(i,j)/wavec(i,j)
        END DO
      END DO
!
# elif defined ROLLER_RENIERS
!
!  Solve roller evolution equation for rollA.
!
      numits=30
      DO it=1,numits
!
!  Compute roller breaking source term (Eqn 40) and 
!  roller disspation sink term (Eqn 41).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff3=(ubar(i,j,nrhs)+wavenx(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*on_u(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FX(i,j)=cff1*rollA(i-1,j)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff3=(vbar(i,j,nrhs)+waveny(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*om_v(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FE(i,j)=cff1*rollA(i,j-1)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)/REAL(numits,r8)
            cff1=cff*(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
            rollA(i,j)=rollA(i,j)-cff1
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Dissip_roller(i,j)=g*sinb*rollA(i,j)*sigma(i,j)/wavec(i,j) 
          END DO
        END DO
!
!  Add roller source / sink term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)/REAL(numits,r8)
            rollA(i,j)=rollA(i,j)+cff*osigma(i,j)*                      &
     &                 (wec_alpha(ng)*Dissip_break(i,j)-                &
     &                 Dissip_roller(i,j))
          END DO
        END DO
# endif
!
!  Call bc's.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    rollA)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      rollA)
# endif
# if defined ROLLER_RENIERS
      END DO
# endif
!
!  Apply boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Dissip_roller)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_roller)
# endif
      RETURN
      END SUBROUTINE wec_roller_tile
#endif
      END MODULE wec_roller_mod
