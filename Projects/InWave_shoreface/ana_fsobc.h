      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id: ana_fsobc.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff, fac, omega, phase, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined INLET_TEST
      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        cff=-1.0_r8*sin(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_north(i)=cff
        END DO
      END IF
#elif defined KELVIN
      fac=1.0_r8                                ! zeta0
      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Istr-1,j)/      &
     &                SQRT(g*GRID(ng)%h(Istr-1,j)))
          BOUNDARY(ng)%zeta_west(j)=val*COS(omega*time(ng))
        END DO
      END IF

      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          cff=1.0_r8/SQRT(g*GRID(ng)%h(Istr-1,j))
          val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Iend,j)*cff)
          BOUNDARY(ng)%zeta_east(j)=val*COS(omega*GRID(ng)%xp(Iend,j)*  &
     &                                      cff-omega*time(ng))
        END DO
      END IF
#elif defined ESTUARY_TEST
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        cff=1.0_r8*SIN(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
#elif defined SED_TEST1
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        fac=100.0_r8
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=9.0E-06_r8*fac
        END DO
      END IF

      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        fac=100.0_r8
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=9.0E-06_r8*REAL(Iend+1,r8)*fac
        END DO
      END IF
#elif defined INWAVE_SHOREFACE
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
!!      cff=-1.0_r8*SIN(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        cff=0.0_r8
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
#elif defined TEST_CHAN
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        cff=0.0_r8
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF

      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        cff=-0.4040_r8*MIN(time(ng)/150000.0_r8,1.0_r8)
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=cff
        END DO
      END IF
#elif defined WEDDELL
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)/REAL(Iend+1,r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8)/REAL(Iend+1,r8))*deg2rad
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=fac*val*COS(omega-phase)
        END DO
      END IF

      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8))*deg2rad
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=fac*val*COS(omega-phase)
        END DO
      END IF
#else
      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF

      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF

      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
