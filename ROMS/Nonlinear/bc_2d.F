#include "cppdefs.h"
      MODULE bc_2d_mod
!
!svn $Id: bc_2d.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 2D fields.                                                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r2d_tile        Boundary conditions for field at RHO-points    !
!    bc_u2d_tile        Boundary conditions for field at U-points      !
!    bc_v2d_tile        Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
! 
!***********************************************************************
      SUBROUTINE bc_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

#if !defined EW_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Iend+1,j)=A(Iend,j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Istr-1,j)=A(Istr,j)
        END DO
      END IF
#endif

#if !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Jend+1)=A(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Jstr-1)=A(i,Jstr)
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                      &
     &                           A(Istr-1,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE bc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

#if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

#if !defined EW_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient 
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
# ifdef EASTERN_WALL
          A(Iend+1,j)=0.0_r8
# else
          A(Iend+1,j)=A(Iend,j)
# endif
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
# ifdef WESTERN_WALL
          A(Istr,j)=0.0_r8
# else
          A(Istr,j)=A(Istr+1,j)
# endif
        END DO
      END IF
#endif

#if !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
# ifdef NORTHERN_WALL
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          A(i,Jend+1)=gamma2(ng)*A(i,Jend)
#  ifdef MASKING
          A(i,Jend+1)=A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
#  undef I_RANGE
# else
        DO i=IstrU,Iend
          A(i,Jend+1)=A(i,Jend)
        END DO
# endif
      END IF

      IF (SOUTHERN_EDGE) THEN
# ifdef SOUTHERN_WALL
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          A(i,Jstr-1)=gamma2(ng)*A(i,Jstr)
#  ifdef MASKING
          A(i,Jstr-1)=A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
#  undef I_RANGE
# else
        DO i=IstrU,Iend
          A(i,Jstr-1)=A(i,Jstr)
        END DO
# endif
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                      &
     &                           A(Istr  ,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                      &
     &                           A(Istr+1,Jend+1)) 
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE bc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

#if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

#if !defined EW_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
# ifdef EASTERN_WALL
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          A(Iend+1,j)=gamma2(ng)*A(Iend,j)
#  ifdef MASKING
          A(Iend+1,j)=A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO
#  undef J_RANGE
# else
        DO j=JstrV,Jend
          A(Iend+1,j)=A(Iend,j)
        END DO
# endif
      END IF

      IF (WESTERN_EDGE) THEN
# ifdef WESTERN_WALL
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          A(Istr-1,j)=gamma2(ng)*A(Istr,j)
#  ifdef MASKING
          A(Istr-1,j)=A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#  endif
        END DO
#  undef J_RANGE
# else
        DO j=JstrV,Jend
          A(Istr-1,j)=A(Istr,j)
        END DO
# endif
      END IF
#endif

#if !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
# ifdef NORTHERN_WALL
          A(i,Jend+1)=0.0_r8
# else
          A(i,Jend+1)=A(i,Jend)
# endif
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
# ifdef SOUTHERN_WALL
          A(i,Jstr)=0.0_r8
# else
          A(i,Jstr)=A(i,Jstr+1)
# endif
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                      &
     &                           A(Istr-1,Jstr+1))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                      &
     &                           A(Iend+1,Jstr+1))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE bc_v2d_tile

      END MODULE bc_2d_mod
