      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id: ana_fsobc.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
# ifdef NORTH_FSOBC
      IF (NORTHERN_EDGE) THEN
        cff=-1.0_r8*sin(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_north(i)=cff
        END DO
      END IF
# endif
#elif defined KELVIN
      fac=1.0_r8                                ! zeta0
      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Istr-1,j)/      &
     &                SQRT(g*GRID(ng)%h(Istr-1,j)))
          BOUNDARY(ng)%zeta_west(j)=val*COS(omega*time(ng))
        END DO
      END IF
# endif
# ifdef EAST_FSOBC
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=1.0_r8/SQRT(g*GRID(ng)%h(Istr-1,j))
          val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Iend,j)*cff)
          BOUNDARY(ng)%zeta_east(j)=val*COS(omega*GRID(ng)%xp(Iend,j)*  &
     &                                      cff-omega*time(ng))
        END DO
      END IF
# endif
#elif defined ESTUARY_TEST
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        cff=1.0_r8*SIN(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
# endif
#elif defined SED_TEST1
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        fac=100.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=9.0E-06_r8*fac
        END DO
      END IF
# endif
# ifdef EAST_FSOBC
      IF (EASTERN_EDGE) THEN
        fac=100.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=9.0E-06_r8*REAL(Iend+1,r8)*fac
        END DO
      END IF
# endif
#elif defined SHOREFACE
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
!!      cff=-1.0_r8*SIN(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
        cff=0.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
# endif

#elif defined DUCKNC
#ifdef EAST_FSOBC
     IF (EASTERN_EDGE) THEN
     cff=0.70_r8
      DO j=JstrR,JendR
!!     BOUNDARY(ng)%zeta_east(j)=cff*tanh(time(ng)/500.0_r8)
       BOUNDARY(ng)%zeta_east(j)=cff
      END DO
     END IF
#endif
#elif defined TEST_CHAN
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        cff=0.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=cff
        END DO
      END IF
# endif
# ifdef EAST_FSOBC
      IF (EASTERN_EDGE) THEN
        cff=-0.4040_r8*MIN(time(ng)/150000.0_r8,1.0_r8)
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=cff
        END DO
      END IF
# endif
#elif defined WEDDELL
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)/REAL(Iend+1,r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8)/REAL(Iend+1,r8))*deg2rad
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=fac*val*COS(omega-phase)
        END DO
      END IF
# endif
# ifdef EAST_FSOBC
      IF (EASTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8))*deg2rad
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=fac*val*COS(omega-phase)
        END DO
      END IF
# endif
#else
# ifdef EAST_FSOBC
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF
# endif
# ifdef WEST_FSOBC
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF
# endif
# ifdef SOUTH_FSOBC
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
# endif
# ifdef NORTH_FSOBC
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
# endif
#endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
