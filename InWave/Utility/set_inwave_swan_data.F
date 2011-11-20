#include "cppdefs.h"
#if defined INWAVE_COUPLING & defined INWAVE_SWAN_COUPLING
      SUBROUTINE set_inwave_swan_data (ng, tile)
!
!svn $Id: set_inwave_swan_data.F 799 2009-12-08 20:38:55Z jcwarner $
!=======================================================================
!                                                                      !
!  This routine computes the action density for the boundary condition !
!  The action density is derived from the directional                  !
!  wave spectra given by SWAN                                          !
!=======================================================================
!
      USE mod_param
      USE mod_inwave_bound
!
!  Imported variable declarations.

      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_inwave_swan_data_tile (ng, tile,                         &
     &                           LBi, UBi, LBj, UBj)
# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_inwave_swan_data
!
!***********************************************************************
      SUBROUTINE set_inwave_swan_data_tile (ng, tile,                   &
     &                                 LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_inwave_bound
      USE mod_inwave_swan
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d
!
!  Local variable declarations.
!

      real(r8)::twopi

# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      logical :: update = .FALSE.

# include "set_bounds.h"

# ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif

#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif

      twopi=2.0_r8*pi

# if defined WEST_AC_OBC
      DO d=1,ND(ng)
        DO j=Jmin,Jmax
          WAVEB(ng)%AC_west(j,d)=WAVES(ng)%Ampzeta(iic(ng),d)/          &
     &    WAVEP(ng)% Tr(LBi,j,ND(ng))*twopi
        END DO
      END DO
# endif
# if defined EAST_AC_OBC
      DO d=1,ND(ng)
        DO j=Jmin,Jmax
          WAVEB(ng)%AC_east(j,d)=WAVES(ng)%Ampzeta(iic(ng),d)/          &
     &    WAVEP(ng)% Tr(UBi,j,ND(ng))*twopi
        END DO
      END DO
# endif
# if defined NORTH_AC_OBC
      DO d=1,ND(ng)
        DO i=Imin,Imax
          WAVEB(ng)%AC_north(i,d)=WAVES(ng)%Ampzeta(iic(ng),d)/         &
     &    WAVEP(ng)% Tr(i,UBj,ND(ng))*twopi
        END DO
      END DO
# endif
# if defined SOUTH_AC_OBC
      DO d=1,ND(ng)
        DO j=Imin,Imax
          WAVEB(ng)%AC_south(j,d)=WAVES(ng)%Ampzeta(iic(ng),d)/         &
     &    WAVEP(ng)% Tr(i,LBj,ND(ng))*twopi
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE set_inwave_swan_data_tile
#else
      SUBROUTINE set_inwave_swan_data
      RETURN
      END SUBROUTINE set_inwave_swan_data
#endif
