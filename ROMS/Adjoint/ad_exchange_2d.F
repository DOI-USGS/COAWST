#include "cppdefs.h"
      MODULE ad_exchange_2d_mod
#if defined ADJOINT && (defined EW_PERIODIC || defined NS_PERIODIC)
!
!svn $Id: ad_exchange_2d.F 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package contains adjoint periodic boundary conditions routines !
!  for 2D primitive variables.                                         !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!   ad_exchange_p2d_tile    periodic conditions at PSI-points          !
!   ad_exchange_r2d_tile    periodic conditions at RHO-points          !
!   ad_exchange_u2d_tile    periodic conditions at U-points            !
!   ad_exchange_v2d_tile    periodic conditions at V-points            !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_exchange_p2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A(1 ,Mm(ng)-2)
!>
          ad_A(1 ,Mm(ng)-2)=ad_A(1 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A(1 ,Mm(ng)-1)
!>
          ad_A(1 ,Mm(ng)-1)=ad_A(1 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A(1 ,Mm(ng)  )
!>
          ad_A(1 ,Mm(ng)  )=ad_A(1 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A(2 ,Mm(ng)-2)
!>
          ad_A(2 ,Mm(ng)-2)=ad_A(2 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A(2 ,Mm(ng)-1)
!>
          ad_A(2 ,Mm(ng)-1)=ad_A(2 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A(2 ,Mm(ng)  )
!>
          ad_A(2 ,Mm(ng)  )=ad_A(2 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1, 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A( 2, 1)
!>
          ad_A( 2, 1)=ad_A( 2, 1)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        ENDIF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_r2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A(1  ,1  )
!>
          ad_A(1  ,1  )=ad_A(1  ,1  )+                                  &
     &                  ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A(1  ,2)
!>
          ad_A(1  ,2)=ad_A(1  ,2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A(2,1  )
!>
          ad_A(2,1  )=ad_A(2,1  )+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
!>
          ad_A(2,2)=ad_A(2,2)+                                          &
     &              ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_u2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1, 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A( 2, 1)
!>
          ad_A( 2, 1)=ad_A( 2, 1)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_v2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1 , 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A(2,1  )
!>
          ad_A( 2, 1)=ad_A( 2, 1 )+                                     &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF (NtileJ(ng).eq.1) THEN
# endif
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef I_RANGE
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1  ,j)
!>
            ad_A(1  ,j)=ad_A(1  ,j)+                                    &
     &                  ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF
#  endif
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_v2d_tile
#endif
      END MODULE ad_exchange_2d_mod
