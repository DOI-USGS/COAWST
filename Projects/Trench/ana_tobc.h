      SUBROUTINE ana_tobc (ng, tile, model)
!
!! svn $Id: ana_tobc.h 1328 2008-01-23 03:20:41Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets tracer-type variables open boundary conditions    !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
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
      CALL ana_tobc_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, nstp(ng),                 &
     &                    GRID(ng) % z_r,                               &
     &                    OCEAN(ng) % t)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(34)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tobc
!
!***********************************************************************
      SUBROUTINE ana_tobc_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj, nstp,               &
     &                          z_r, t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_boundary
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, nstp

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, ised, itrc, j, k
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tracers open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined TRENCH
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_west(j,k,itemp)=10.0_r8
!            BOUNDARY(ng)%t_west(j,k,isalt)=0.0_r8
#  if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_west(j,k,idsed(ised))=                     &
     &                                  t(1,j,k,nstp,idsed(ised))
            END DO
#  endif
          END DO
        END DO
      END IF
#else
      IF (EASTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_east(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_west(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_south(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_north(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_tobc_tile
