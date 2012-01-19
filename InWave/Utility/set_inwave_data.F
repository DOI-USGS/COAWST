#include "cppdefs.h"
#ifdef INWAVE_MODEL
      SUBROUTINE set_inwave_data (ng, tile)
!
!svn $Id: set_data.F 799 2008-10-20 20:38:55Z jcwarner $
!
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary input data.             !
!  It time-interpolates between snapshots.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_inwave_data_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_inwave_data
!
!***********************************************************************
      SUBROUTINE set_inwave_data_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_bound
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

# ifdef AC_OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: ND_BND

!
!  Local variable declarations.
!
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

# ifdef AC_OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif

      ND_BND=WAVEB(ng)%ND_BND


# if defined WEST_AC_OBC
        CALL set_ngfld (ng, iNLM, idACbc(iwest), JLB, JUB, ND_BND,      &
     &                  0, Mm(ng)+1, ND_BND,                            &
     &                  WAVEB(ng) % ACG_west(JLB,1,1),                  &
     &                  WAVEB(ng) % AC_west(JLB,1),                     &
     &                  update)
# endif
# if defined EAST_AC_OBC
        CALL set_ngfld (ng, iNLM, idACbc(ieast), JLB, JUB, ND_BND,      &
     &                  0, Mm(ng)+1, ND_BND,                            &
     &                  WAVEB(ng) % ACG_east(JLB,1,1),                  &
     &                  WAVEB(ng) % AC_east(JLB,1),                     &
     &                  update)
# endif
# if defined NORTH_AC_OBC
        CALL set_ngfld (ng, iNLM, idACbc(inorth), ILB, IUB, ND_BND,     &
     &                  0, Lm(ng)+1, ND_BND,                            &
     &                  WAVEB(ng) % ACG_north(ILB,1,1),                 &
     &                  WAVEB(ng) % AC_north(ILB,1),                    &
     &                  update)
# endif
# if defined SOUTH_AC_OBC
        CALL set_ngfld (ng, iNLM, idACbc(isouth), ILB, IUB, ND_BND,     &
     &                  0, Lm(ng)+1, ND_BND,                            &
     &                  WAVEB(ng) % ACG_south(ILB,1,1),                 &
     &                  WAVEB(ng) % AC_south(ILB,1),                    &
     &                  update)
# endif
      RETURN
      END SUBROUTINE set_inwave_data_tile
#else
      SUBROUTINE set_inwave_data
      RETURN
      END SUBROUTINE set_inwave_data
#endif
