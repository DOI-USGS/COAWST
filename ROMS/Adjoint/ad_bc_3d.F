#include "cppdefs.h"
      MODULE ad_bc_3d_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_bc_3d.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 3D adjoint fields.                                          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_bc_r3d_tile     Boundary conditions for field at RHO-points    !
!    ad_bc_u3d_tile     Boundary conditions for field at U-points      !
!    ad_bc_v3d_tile     Boundary conditions for field at V-points      !
!    ad_bc_w3d_tile     Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bc_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_r3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend  ,Jstr-1,k)+          &
!>   &                                  tl_A(Iend+1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jstr-1,k)=0.5_r8*(tl_A(Istr  ,Jstr-1,k)+          &
!>   &                                  tl_A(Istr-1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr-1,k)
          ad_A(Istr  ,Jstr-1,k)=ad_A(Istr  ,Jstr-1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=ad_A(Istr-1,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_r3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_u3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_u3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr  ,Jend+1,k)=0.5_r8*(tl_A(Istr  ,Jend  ,k)+          &
!>   &                                  tl_A(Istr+1,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr  ,Jend+1,k)
          ad_A(Istr  ,Jend  ,k)=ad_A(Istr  ,Jend  ,k)+adfac
          ad_A(Istr+1,Jend+1,k)=ad_A(Istr+1,Jend+1,k)+adfac
          ad_A(Istr  ,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend  ,Jstr-1,k)+          &
!>   &                                  tl_A(Iend+1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr  ,Jstr-1,k)=0.5_r8*(tl_A(Istr+1,Jstr-1,k)+          &
!>   &                                  tl_A(Istr  ,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr  ,Jstr-1,k)
          ad_A(Istr+1,Jstr-1,k)=ad_A(Istr+1,Jstr-1,k)+adfac
          ad_A(Istr  ,Jstr  ,k)=ad_A(Istr  ,Jstr  ,k)+adfac
          ad_A(Istr  ,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed (free-slip/no-slip)
!  or gradient.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
#  ifdef SOUTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=LBk,UBk
          DO i=I_RANGE
#   ifdef MASKING
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr-1,k)*GRID(ng)%umask(i,Jstr-1)
!>
            ad_A(i,Jstr-1,k)=ad_A(i,Jstr-1,k)*GRID(ng)%umask(i,Jstr-1)
#   endif
!>          tl_A(i,Jstr-1,k)=gamma2(ng)*tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+gamma2(ng)*ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
#   undef I_RANGE
#  else
        DO k=LBk,UBk
          DO i=IstrU,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
#  endif
      END IF

      IF (NORTHERN_EDGE) THEN
#  ifdef NORTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=LBk,UBk
          DO i=I_RANGE
#   ifdef MASKING
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend+1,k)*GRID(ng)%umask(i,Jend+1)
!>
            ad_A(i,Jend+1,k)=ad_A(i,Jend+1,k)*GRID(ng)%umask(i,Jend+1)
#   endif
!>          tl_A(i,Jend+1,k)=gamma2(ng)*tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+gamma2(ng)*ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
#   undef I_RANGE
#  else
        DO k=LBk,UBk
          DO i=IstrU,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
#  ifdef WESTERN_WALL
!>          tl_A(Istr,j,k)=0.0_r8
!>
            ad_A(Istr,j,k)=0.0_r8
#  else
!>          tl_A(Istr,j,k)=tl_A(Istr+1,j,k)
!>
            ad_A(Istr+1,j,k)=ad_A(Istr+1,j,k)+ad_A(Istr,j,k)
            ad_A(Istr  ,j,k)=0.0_r8
#  endif
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
#  ifdef EASTERN_WALL
!>          tl_A(Iend+1,j,k)=0.0_r8
!>
            ad_A(Iend+1,j,k)=0.0_r8
#  else
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
#  endif
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_u3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_v3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_v3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jstr  ,k)=0.5_r8*(tl_A(Iend  ,Jstr  ,k)+          &
!>   &                                  tl_A(Iend+1,Jstr+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr  ,k)
          ad_A(Iend  ,Jstr  ,k)=ad_A(Iend  ,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr+1,k)=ad_A(Iend+1,Jstr+1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jstr  ,k)=0.5_r8*(tl_A(Istr  ,Jstr  ,k)+          &
!>   &                                  tl_A(Istr-1,Jstr+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr  ,k)
          ad_A(Istr  ,Jstr  ,k)=ad_A(Istr  ,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr+1,k)=ad_A(Istr-1,Jstr+1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
#  ifdef SOUTHERN_WALL
!>          tl_A(i,Jstr,k)=0.0_r8
!>
            ad_A(i,Jstr,k)=0.0_r8
#  else
!>          tl_A(i,Jstr,k)=tl_A(i,Jstr+1,k)
!>
            ad_A(i,Jstr+1,k)=ad_A(i,Jstr+1,k)+ad_A(i,Jstr,k)
            ad_A(i,Jstr  ,k)=0.0_r8
#  endif
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
#  ifdef NORTHERN_WALL
!>          tl_A(i,Jend+1,k)=0.0_r8
!>
            ad_A(i,Jend+1,k)=0.0_r8
#  else
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
#  endif
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
#  ifdef WESTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=LBk,UBk
          DO j=J_RANGE
#   ifdef MASKING
!>          tl_A(Istr-1,j,k)=tl_A(Istr-1,j,k)*GRID(ng)%vmask(Istr-1,j)
!>
            ad_A(Istr-1,j,k)=ad_A(Istr-1,j,k)*GRID(ng)%vmask(Istr-1,j)
#   endif
!>          tl_A(Istr-1,j,k)=gamma2(ng)*tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+gamma2(ng)*ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
#   undef J_RANGE
#  else
        DO k=LBk,UBk
          DO j=JstrV,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
#  endif
      END IF

      IF (EASTERN_EDGE) THEN
#  ifdef EASTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=LBk,UBk
          DO j=J_RANGE
#   ifdef MASKING
!>          tl_A(Iend+1,j,k)=tl_A(Iend+1,j,k)*GRID(ng)%vmask(Iend+1,j)
!>
            ad_A(Iend+1,j,k)=ad_A(Iend+1,j,k)*GRID(ng)%vmask(Iend+1,j)
#   endif
!>          tl_A(Iend+1,j,k)=gamma2(ng)*tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+gamma2(ng)*ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
#   undef J_RANGE
#  else
        DO k=LBk,UBk
          DO j=JstrV,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_v3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_w3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_w3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend+1,Jstr  ,k)+          &
!>                                      tl_A(Iend  ,Jstr-1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=LBk,UBk
!>        tl_A(Istr-1,Jstr-1,k)=0.5_r8*(tl_A(Istr  ,Jstr-1,k)+          &
!>   &                                  tl_A(Istr-1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr-1,k)
          ad_A(Istr  ,Jstr-1,k)=ad_A(Istr  ,Jstr-1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=ad_A(Istr-1,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ad_bc_w3d_tile

#endif
      END MODULE ad_bc_3d_mod
