      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id: ana_m2obc.h 735 2007-04-27 14:00:46Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
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
        WRITE (ANANAME(12),'(a,a)') TRIM(Adir), '/ana_m2obc.h'
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

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_east(j)=???
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=???
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_west(j)=???
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_west(j)=???
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_south(i)=???
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_south(i)=???
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_north(i)=???
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_north(i)=???
        END DO
      END IF
#elif defined SANDWAVE
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=0.5_r8*(zeta(Istr-1,j,knew)+h(Istr-1,j)+                  &
     &                zeta(Istr  ,j,knew)+h(Istr  ,j))
          BOUNDARY(ng)%ubar_west(j)=-10.0_r8/val
        END DO
        DO j=Jstr,JendR
          vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                  &
     &                zeta(Iend+1,j,knew)+h(Iend+1,j))
          BOUNDARY(ng)%ubar_east(j)=-10.0_r8/val
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
#else
      ana_m2obc_user.h: No values provided for BOUNDARY(ng)%ubar_xxxx
                                               BOUNDARY(ng)%vbar_xxxx
#endif

      RETURN
      END SUBROUTINE ana_m2obc_tile
