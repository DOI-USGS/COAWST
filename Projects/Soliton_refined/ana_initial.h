      SUBROUTINE ana_initial (ng, tile, model)
!
!! svn $Id: ana_initial.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets initial conditions for momentum and tracer     !
!  type variables using analytical expressions.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      IF (model.eq.iNLM) THEN
        CALL ana_NLMinitial_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            GRID(ng) % h,                         &
#ifdef SPHERICAL
     &                            GRID(ng) % lonr,                      &
     &                            GRID(ng) % latr,                      &
#else
     &                            GRID(ng) % xr,                        &
     &                            GRID(ng) % yr,                        &
#endif
#ifdef SOLVE3D
     &                            GRID(ng) % z_r,                       &
     &                            OCEAN(ng) % u,                        &
     &                            OCEAN(ng) % v,                        &
     &                            OCEAN(ng) % t,                        &
#endif
     &                            OCEAN(ng) % ubar,                     &
     &                            OCEAN(ng) % vbar,                     &
     &                            OCEAN(ng) % zeta)
#ifdef TANGENT
      ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
        CALL ana_TLMinitial_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            kstp(ng),                             &
# ifdef SOLVE3D
     &                            nstp(ng),                             &
     &                            OCEAN(ng) % tl_u,                     &
     &                            OCEAN(ng) % tl_v,                     &
     &                            OCEAN(ng) % tl_t,                     &
# endif
     &                            OCEAN(ng) % tl_ubar,                  &
     &                            OCEAN(ng) % tl_vbar,                  &
     &                            OCEAN(ng) % tl_zeta)
#endif
#ifdef ADJOINT
      ELSE IF (model.eq.iADM) THEN
        CALL ana_ADMinitial_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            knew(ng),                             &
# ifdef SOLVE3D
     &                            nstp(ng),                             &
     &                            OCEAN(ng) % ad_u,                     &
     &                            OCEAN(ng) % ad_v,                     &
     &                            OCEAN(ng) % ad_t,                     &
# endif
     &                            OCEAN(ng) % ad_ubar,                  &
     &                            OCEAN(ng) % ad_vbar,                  &
     &                            OCEAN(ng) % ad_zeta)
#endif
      END IF
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(10)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_initial
!
!***********************************************************************
      SUBROUTINE ana_NLMinitial_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                h,                                &
#ifdef SPHERICAL
     &                                lonr, latr,                       &
#else
     &                                xr, yr,                           &
#endif
#ifdef SOLVE3D
     &                                z_r,                              &
     &                                u, v, t,                          &
#endif
     &                                ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(out) :: u(LBi:,LBj:,:,:)
      real(r8), intent(out) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: t(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
#else
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(out) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
      real(r8), intent(out) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Iless, Iplus, i, itrc, j, k
      real(r8) :: depth, dx, val1, val2, val3, val4, x, x0, y, y0

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initial conditions for 2D momentum (m/s) components.
!-----------------------------------------------------------------------
!
#if defined SOLITON_REFINED
      g=1.0_r8
      x0=40.0_r8
      y0=8.0_r8
      val1=0.395_r8
      val2=0.771_r8*(val1*val1)
      DO j=JstrT,JendT
        DO i=IstrTU,IendT
          x=0.5_r8*(xr(i-1,j)+xr(i,j))-x0
          y=0.5_r8*(yr(i-1,j)+yr(i,j))-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          ubar(i,j,1)=0.25_r8*val4*(6.0_r8*y*y-9.0_r8)*                 &
     &                EXP(-0.5_r8*y*y)
        END DO
      END DO
      DO j=JstrTV,JendT
        DO i=IstrT,IendT
          x=0.5_r8*(xr(i,j-1)+xr(i,j))-x0
          y=0.5_r8*(yr(i,j-1)+yr(i,j))-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          vbar(i,j,1)=2.0_r8*val4*y*(-2.0_r8*val1*TANH(val1*x))*        &
     &                EXP(-0.5_r8*y*y)
        END DO
      END DO
#else
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,1)=0.0_r8
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,1)=0.0_r8
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initial conditions for free-surface (m).
!-----------------------------------------------------------------------
!
#if defined SOLITON_REFINED
      x0=40.0_r8
      y0=8.0_r8
      val1=0.395_r8
      val2=0.771_r8*(val1*val1)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          x=xr(i,j)-x0
          y=yr(i,j)-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          zeta(i,j,1)=0.25_r8*val4*(6.0_r8*y*y+3.0_r8)*                 &
     &                EXP(-0.5_r8*y*y)
        END DO
      END DO
#else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=0.0_r8
        END DO
      END DO
#endif
#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initial conditions for 3D momentum components (m/s).
!-----------------------------------------------------------------------
!
# if defined RIVERPLUME2
      DO k=1,N(ng)
       DO j=JstrR,JendR
         DO i=Istr,IendR
            u(i,j,k,1)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,1)=-0.05_r8*LOG((h(i,j)+z_r(i,j,k))/Zob(ng))/        &
     &                 (LOG(h(i,j)/Zob(ng))-1.0_r8+Zob(ng)/h(i,j))
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)
       DO j=JstrR,JendR
         DO i=Istr,IendR
            u(i,j,k,1)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,1)=0.0_r8
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Initial conditions for tracer type variables.
!-----------------------------------------------------------------------
!
!  Set initial conditions for potential temperature (Celsius) and
!  salinity (PSU).
!
# if defined BENCHMARK
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*800.0_r8/g)*(5.0E-05_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=val2*EXP(z_r(i,j,k)/800.0_r8)*             &
     &                       (0.6_r8-0.4_r8*TANH(z_r(i,j,k)/800.0_r8))
            t(i,j,k,1,isalt)=35.0_r8
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)
#  ifdef SALINITY
            t(i,j,k,1,isalt)=S0(ng)
#  endif
          END DO
        END DO
      END DO
# endif
#endif
      RETURN
      END SUBROUTINE ana_NLMinitial_tile

