#include "cppdefs.h"
#ifdef INWAVE_COUPLING
      SUBROUTINE set_inwave_data (ng, tile)
!
!svn $Id: set_data.F 799 2008-10-20 20:38:55Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  assimilation input data. It time-interpolates between snapshots.    !
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
      USE mod_boundary
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
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
# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
#endif
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

# ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif

!      IF (SOUTH_WEST_TEST) THEN
!#    ifdef WEST_M3OBC


        CALL set_ngfld (ng, iNLM, idACbc(iwest), JLB, JUB, ND(ng),      &
     &                  0, Mm(ng)+1, ND(ng),                            &
     &                  WAVEB(ng) % ACG_west(JLB,1,1),                  &
     &                  WAVEB(ng) % AC_west(JLB,1),                     &
     &                  update)

        write(*,*) WAVEB(ng) % AC_west(1,1)


      RETURN
      END SUBROUTINE set_inwave_data_tile
#else
      SUBROUTINE set_inwave_data
      RETURN
      END SUBROUTINE set_inwave_data
#endif
