      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id: ana_m2obc.h 1328 2008-01-23 03:20:41Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
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
      CALL ana_m2obc_tile (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(12)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           knew,                                  &
     &                           angler, h, pm, pn, on_u,               &
#ifdef MASKING
     &                           umask,                                 &
#endif
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)   
      real(r8), intent(in) :: h(LBi:,LBj:)   
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)   
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)   
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: angle, cff, fac, major, minor, omega, phase, val
      real(r8) :: ramp
#if defined TRENCH
      real(r8) :: my_area, my_width
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined TRENCH
      IF (WESTERN_EDGE) THEN
        my_area=0.0_r8
        my_width=0.0_r8
        DO j=Jstr,Jend
          my_area=my_area+0.5_r8*(zeta(Istr-1,j,knew)+h(Istr-1,j)+      &
     &                            zeta(Istr  ,j,knew)+h(Istr  ,j))*     &
     &                           on_u(Istr,j)
          my_width=my_width+on_u(Istr,j)
        END DO
        fac=my_width*0.39_r8*0.51_r8               !(width  depth  ubar)
        DO j=Jstr,Jend
          BOUNDARY(ng)%ubar_west(j)=fac/my_area
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        my_area=0.0_r8
        my_width=0.0_r8
        DO j=Jstr,Jend
          my_area=my_area+0.5_r8*(zeta(Iend+1,j,knew)+h(Iend+1,j)+      &
     &                            zeta(Iend  ,j,knew)+h(Iend  ,j))*     &
     &                           on_u(Iend,j)
         my_width=my_width+on_u(Iend,j)
        END DO
        fac=my_width*0.39_r8*0.51_r8               !(width  depth  ubar)
        DO j=Jstr,Jend
          BOUNDARY(ng)%ubar_east(j)=fac/my_area
        END DO
      END IF
#else
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_west(j)=0.0_r8
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_south(i)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_m2obc_tile
